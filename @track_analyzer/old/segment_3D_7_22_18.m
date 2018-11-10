% 3D Segmentation of cell nuclei. Based on: 'Object segmentation and ground
% truth in 3D embryonic imaging'. 

% Specify segmentation parameters in params, and function options in step.
% Written for single time-point images. <<< will write for looping over
% time later...

function obj = segment_3D(obj, params, step)



%% Pre-processing using bio-formats. 

%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = bfGetReader(obj.exp_info.img_file);
series = 1;


%Get the image size of this series. 

T = reader.getSizeT;

 %Get step.debug
if(step.debug)
    params
    T = 1;
end  

%% Generate im_info structure
    ZSlicesinStack = obj.exp_info.z_planes;

    image_bits     = reader.getBitsPerPixel;   
    
    %Since we do parallel processing for time points, pre-compute the list
    %of plane indices for each time point. 
    Z_cell = {};
    for t = 1:T
        planes = [];
        for Z = 1:ZSlicesinStack
            planes(Z) = reader.getIndex(Z-1,params.seg_channel-1,t-1)+1;
        end
        Z_cell{t} = planes;
    end
    %Close
    reader.close();
    
    
%% Loop over images. Need to check if more than one image, do parallel. 
if(T>1)
    parallel = 0;
else
    parallel = 0;
end

%Experiment info to pass along. 
exp_info = obj.exp_info;

%Now check to see if any frames already exist. 
frame_files = obj.get_frame_files;
if(isempty(frame_files{1}) || step.FORCE_ALL_FRAMES)
    new_ts = [1:T];
else
    for i = 1:length(frame_files)
        [~,fname,~] = fileparts(frame_files{i});
        t_val(i) = str2num(fname(end-3:end));
    end
    all_ts = [1:T];
    new_ts = setdiff(all_ts,t_val);
end


%Now run loops 
    if(parallel)
        disp('New frames to calculate:')
        disp(new_ts)
        parfor g = 1:length(new_ts)
            t = new_ts(g);
            inner_function( exp_info, Z_cell{t}, params, step, t)
        end
    else
        disp('New frames to calculate:')
        disp(num2str(new_ts))
        for t = new_ts
            inner_function( exp_info, Z_cell{t},params,step,t)
        end
    end


end

%% Inner function to run 
function inner_function( exp_info, planes, params, step, t)
    
%Display 
disp(['Started frame: ',num2str(t)])

%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%  
    %% Segmentation steps
        display_plots = step.display_plots;


    %% Set Segmentation Algorithm Parameter values
    
        AbsMinVol = params.AbsMinVol;
        AbsMaxVol = params.AbsMaxVol; %Maximum voxel per single nucleus

        %%Derivatives Sum (DS) Algorithm
        alpha = params.alpha;
        beta = params.beta;
        gamma = params.gamma;
        epsilon =params.epsilon;
        delta = params.delta;
        sigmagradient=params.sigmagradient;

        %%Noise filter parameters
        h = params.noise_filter; 
        w1=params.median_filter(1);w2=params.median_filter(2); %%Median filter

        %%Perona-Malik based non-linear isotropic diffusion filter
        diffuse_iterations=params.diffuse_iterations;
        kappa1=params.kappa1;
        kappa2=params.kappa2;
        option=params.option;

        %Thresholding
        thrshlevel = params.thrshlevel;
        % nconn value for neighborhood calculations
        nconn_BW = params.nconn_BW;
        nconn_BW2= params.nconn_BW2;

        MeanFilterSensitivity = params.MeanFilterSensitivity;
        MeanFilterNeighborhood = params.MeanFilterNeighborhood;

        xscale         = exp_info.pixel_size(1);
        yscale         = exp_info.pixel_size(2);
        zscale         = exp_info.pixel_size(3);

%%%%%%%%%%%%%%% START PROCESSING %%%%%%%%%%%%%%%%
    reader = bfGetReader(exp_info.img_file);
    image_bits     = reader.getBitsPerPixel;
    size_x = reader.getSizeX;
    size_y = reader.getSizeY;
    ZSlicesinStack = length(planes);

    %Get image planes for this time point. 
    I = zeros(size_y,size_x,ZSlicesinStack);

    %Get the bio-formats image index corresponding to this z-stack:
    for i = 1:ZSlicesinStack
        this_plane_img = bfGetPlane(reader,planes(i));
        I(:,:,i)     = this_plane_img;
    end
    %Close reader
    reader.close();
    
    %Reduce i for debugging
    if(step.debug)
        %I = I(1:300,1:300,:);
    end
    Ie = zeros(size(I));
    G = Ie; Diff_im = Ie; Fim = Ie; Mag_fim = Ie; Mag_Lim = Ie; Det_hessian = Ie;
    Multi = Ie; J = Ie;
     
    
    %% PRE-PROCESSING STEPS. Filters, and removing outliers. 
       
    %PRE-PROCESS image with mean filter if desired
    if(step.MeanFilter)
        tmp = adaptthresh(uint16(I),MeanFilterSensitivity,'NeighborhoodSize',MeanFilterNeighborhood);
        
        %%%% This step might help bring everything to 16bit levels. Try
        %%%% manually changing bit size to 16 here. 
        I = tmp.*65000;
        image_bits = 16;
    end
    
    G = medfilt3(I, [w1 w2 w3]);

    
    %% Diffusion 
     for zslice= 1:ZSlicesinStack %% iterate through number of slices
        %
        %use imshow to plot any of the Steps from 1-5. eg.figure; imshow((I(:,:,zslice))

        %%Step1: Denoising Filters
        %G(:,:,zslice) = imfilter(I(:,:,zslice), h,'replicate');
        %%Optional filters median or deconvlucy
        %G(:,:,zslice) = deconvlucy(G(:,:,zslice),h); %%use the same 'h' as for Gaussian or define new variable

        %%Step2: Perona & Malik non-linear isotropic diffusion filter, refer to
        %%diffusioncode.m for details, here alter only diffuse_iterations & kappa1       
        Diff_im(:,:,zslice) = diffusioncode(I(:,:,zslice), diffuse_iterations, 0.1429, kappa1, kappa2, option);

        %%Step3a: Gauss gradient, 1st derivative of the image, refer to gaussgradient.m for details
        Fim(:,:,zslice)=mat2gray(Diff_im(:,:,zslice));
        [imx,imy]=gaussgradient(Fim(:,:,zslice),sigmagradient);
        Mag_fim(:,:,zslice) = hypot(imx,imy);

        %%Step3b: Laplacian
        [L_imxx, L_imxy] = gaussgradient(imx,sigmagradient);
        [L_imyx, L_imyy] = gaussgradient(imy,sigmagradient);
        Mag_Lim(:,:,zslice) = L_imxx+L_imyy;  %Laplacian (trace of the 2D matrix)
        Mag_Lim(:,:,zslice)=Mag_Lim(:,:,zslice).*(Mag_Lim(:,:,zslice)>0);

        %%Step3c: Hessian Determniant
        Det_hessian(:,:,zslice)=L_imxx.*L_imyy-L_imxy.*L_imyx;
        Det_hessian(:,:,zslice)=Det_hessian(:,:,zslice).*(Det_hessian(:,:,zslice)<0);
        %%Step4: Masking function using tanh on the Summed Derivatives from Step3
        X=gamma-((alpha*Mag_fim(:,:,zslice) + beta*Mag_Lim(:,:,zslice) + epsilon*abs(Det_hessian(:,:,zslice))) /delta);
        Multi(:,:,zslice)=0.5*(tanh(X)+1);% masking function
        J(:,:,zslice)=(double(G(:,:,zslice)).*Multi(:,:,zslice)); % masked image applied on smoothened image
    end

        if(step.debug)
            imshow3D(J)
            pause
        end

        %% Global Otsu's method. Image Thresholding using Otsu's method, 8-bit or 16-bit

        BW = zeros(size(I));
        if (image_bits==8)
            %thrshlevel = otsuthresh( uint8(J(:)));
            BW(:,:,zslice)=im2bw(uint8(J(:,:,zslice)),thrshlevel);
        elseif(image_bits==16)
            %thrshlevel = otsuthresh( uint16(J(:)));
            for i = 1:1:ZSlicesinStack
                %Threshold 
                BW(:,:,i) = im2bw(uint16(J(:,:,i)),thrshlevel);
                %Also fill holes laterally. Vertical holes can be
                %problematic if we just to imfill with 3D stack
                BW(:,:,i) = imfill(BW(:,:,i),'holes');
            end

        end
        
    %% Looking for high intensity outliers. 
    if(step.FilterOutliers)
       
        %Estimate background value. 
        bg_img = I(~BW);
        mean_bg_value = mean(bg_img(:));
        
        %Find outliers using a high threshold. 
        high = I > params.OutlierThreshold;
        
        %Dilate 
        se = strel('cuboid',params.MeanFilterNeighborhood );
        dil = imdilate(high, se);
        
        %Collect region info. 
        stats = regionprops(logical(dil),'PixelIdxList','Image');
        
        %Change high-intensity regions to NaN or mean bg value? 
        all_idx = cat(1,stats.PixelIdxList);
        I(all_idx) = mean_bg_value;
        
        %Fill in missing data with inpaintn
        %I = inpaintn(I);
        
        %Alternatively, we can fill in with background noise?
        if(step.debug)
            imshow3D(I);
            pause
        end
        
        BW(all_idx) = 0;
    end
    
    %% Step6: Connect similar pixels of zslices to get a 3D binary stack

    %%Step7: Computes three properties for each 3D segmented object
    Con=bwconncomp(logical(BW),nconn_BW);
    stats = regionprops(Con,'Centroid','Area','PixelList','PixelIdxList');
    numobj=numel(stats); %total number of segmented objects
    
    if(step.debug)
        imshow3D(BW)
        pause
    end
    %%Step8: Image stack dimensions & image scaling in mirons
    [nx,ny,nz]=size(I);
    x=(1:nx)*xscale;
    y=(1:ny)*yscale;
    z=(1:nz)*zscale;
        
    %% Figure out which volumes are above the threshold. These might be fused objects. 
    
    %Remove volumes below size threshold. 
    volumes = [stats.Area]';
    sel = volumes < AbsMinVol;
    stats = stats(sel);
    
    %New volume list
    volumes = [stats.Area]';   
    liar_list = find( volumes > AbsMaxVol);
    
    %% Post-Processing of Fused Objects with GMM & Kmeans
     %%By setting to zero, you can switch of post-processing of fused nuclei
    if step.fused_nuclei
        
        %%setting scales to find next peak for fused nuclei based on
        %%image resolution in x,y,z (assessed from micrsocope resolution)
        %%used for voxel frequency distribution
        scale_peak_distx=round(4/xscale);
        scale_peak_disty=round(4/yscale);
        scale_peak_distz=round(3/zscale);
        
        noteindex_GMM =[];%%initialize a variable to track indices of objects that remain fused after GMM,
        %%in case many are fused after GMM, then change parametrs of noise
        %%filters & DS algorithm to improve segmentation accuarcy
        
        %Make a marker mask to create regional minima for seeded 3D
        %watershed. 
        marker = zeros(size(BW));
        
        %%post-process each fused object from liar_list variable
        for l=1:length(liar_list)

            liar_Area=stats(liar_list(l)).Area;
            Pxlist=stats(liar_list(l)).PixelList;
            PxIdx = stats(liar_list(l)).PixelIdxList;
            Px=Pxlist;
            
            %%micron scaling of pixels
            Px(:,1)=Px(:,1)*xscale;
            Px(:,2)=Px(:,2)*yscale;
            Px(:,3)=Px(:,3)*zscale;
            
            %%micron scaling of Centroids
            X_Centroid_px=(stats(liar_list(l)).Centroid(1))*xscale;
            Y_Centroid_px=(stats(liar_list(l)).Centroid(2))*yscale;
            Z_Centroid_px=(stats(liar_list(l)).Centroid(3))*zscale;
            
            Pxlist(end+1,1:3)=[-1 -1 -1]; % to avoid index out of bound error
            
            %%Determine frequency distribution of voxels in x, y, z
            %%directions.Use inter-peak distances within voxel frequency
            %%distribution to determine number of fused nuclei
            
            %%X-Direction frequency distribution
            freqx=tabulate(Pxlist(:,1));
            if (length(freqx(:,1))>10)
                if (length(freqx(2:end,2))>5)
                    peak_distx=scale_peak_distx;
                else
                    peak_distx=1;
                end
                
                if peak_distx==1
                    x_peaks=0;
                else
                    x_peaks=numel(findpeaks(freqx(2:end,2),'minpeakdistance',peak_distx));
                end
            else
                x_peaks=0;
            end
            
            %%Y-Direction frequency distribution
            freqy=tabulate(Pxlist(:,2));
            if (length(freqy(:,1))>10)
                if (length(freqy(2:end,2))>5), peak_disty=scale_peak_disty;else peak_disty=1;end
                if peak_disty==1
                    y_peaks=0;
                else
                    y_peaks=numel(findpeaks(freqy(2:end,2),'minpeakdistance',peak_disty));
                end
            else
                y_peaks=0;
            end
            
            %%Z-Direction frequency distribution
            freqz=tabulate(Pxlist(:,3));
            if (length(freqz(:,1))>5)
                peak_distz=scale_peak_distz;
                if (length(freqz(2:end,2))>5),z_peaks=numel(findpeaks(freqz(2:end,2),'minpeakdistance',peak_distz)); else z_peaks=0;end
            else
                z_peaks=0;
            end
            
            %% Plot VoxelFrequencyDistributionPlot
            if step.display_plots
                VoxelFrequencyDistributionPlot(I,freqx,freqy,freqz); %% Plot of frequency distribution of voxels here
                subplot(3,1,1);title('Voxel frequency distribution for finding peaks for GMM & Kmeans','FontSize',10);
                
            end
            %%
            
            %% Render nuclei in 3D for the respective nuclei under post-processing
            if step.display_plots
                RenderNucleiDS(G,x,y,z,freqx,freqy,freqz,X_Centroid_px,Y_Centroid_px,Z_Centroid_px);
                pause
            end
            
            %%determine possible number of fused nuclei (maximum of peaks obtained from x,y,z directions)
            nuclei_cluster=max([x_peaks y_peaks z_peaks])
            n_cluster(l) = nuclei_cluster;
            
            %% Kmeans Post-Processing/Local Clustering of fused Candidate Nuclei using Kmeans
            %By setting to zero, you can switch off Kmeans
            IDX=0;% for Kmeans to check number of clusters
            if step.Active_Kmeans
                
                if (nuclei_cluster>=1)
                    
                    %Perform K-means clustering. 
                    [IDX,C] = kmeans(Px,nuclei_cluster,'distance','sqeuclidean','onlinephase','on','emptyaction','drop','replicates',10);
                    %Get pixel coordinates of the cluster centers. 
                    C_px = round(C./exp_info.pixel_size);
                    C_px_idx = sub2ind(size(BW),C_px(:,2),C_px(:,1),C_px(:,3));
                    
                    %Add markers to marker binary
                    marker(C_px_idx) = 1;
                end
            end % end Kmeans clsuering
                        
        end
    
        
    end % end of fused nuclei processing

    %% Apply 3D watershedding. 
    if( step.watershed3D)
        tic
        disp('Starting watershed')
        %Make bw of clusters only
        bw_clust = zeros(size(BW));
        idx = cat(1,stats(liar_list).PixelIdxList);
        bw_clust(idx) = 1;
        %Anisotropic distance transform
        imgDist = -bwdistsc(~bw_clust,exp_info.pixel_size);
        %Smoothing
        imgDist = medfilt3(imgDist,[5,5,5]);
        %Get seeds 
        mask = imextendedmin(imgDist,0.7,6); %Seems like smaller neighborhood works better?
        imgDist = imimposemin(imgDist,mask);
        imgDist(~bw_clust) = -inf;
        %Perform marked watershed
        imgLabel = watershed(imgDist);
        %Get rid of '1' label. 
        sel= imgLabel==1;
        imgLabel(sel) = 0;
        stats_fused = regionprops(imgLabel,'Centroid','Area','PixelList','PixelIdxList');
        stats_fused = stats_fused(2:end);
        %Now get rid of original regions that were fused. 
        keep = setdiff([1:length(stats)],liar_list);
        stats = stats(keep);
        %Now append the fused nuclei stats. 
        stats = [stats;stats_fused];
        %Now we need to adjust BW... 
            %Clear out clustered regions
            BW(idx) = 0;
            %Add in new regions
            new_idx = cat(1,stats_fused.PixelIdxList);
            BW(new_idx) = 1;
        disp('Finished watershed')
        toc
    end
    
    if(step.debug)
        figure
       imshow3D(BW) 
    end
    %% Error if there's no cells.
    if(~isfield(stats,'PixelList'))
        error('No cells found!')
    end

        %% Create frame_obj and save. 
        frame_obj = [];
        %Add pixel list/centroid for each region to frame_obj
        counter = 1;
        for i = 1:length(stats)
            if (length(stats(i).PixelIdxList) >= AbsMinVol)
                frame_obj.PixelIdxList{counter} = stats(i).PixelIdxList;
                frame_obj.centroids{counter}    = stats(i).Centroid;
                counter = counter + 1;
            end
        end

        %Add final binarized image to frame_obj for save keeping
        frame_obj.BW = BW;

        %Save frame_obj as done before. 
        fname = ['frame_',sprintf('%04d',t),'.mat'];
        if(~exist(exp_info.nuc_seg_dir,'dir'))
            mkdir(exp_info.nuc_seg_dir)
        end
        parsave([exp_info.nuc_seg_dir,fname],frame_obj)
        disp(['Finished frame:     ',num2str(t)])
end


%Parallel saving technique
function parsave(fname, frame_obj)
save(fname, 'frame_obj')
end

