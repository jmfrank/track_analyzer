% f3D Segmentation of cell nuclei. Based on: 'Object segmentation and ground
% truth in 3D embryonic imaging'. 

% Specify segmentation parameters in params, and function options in step.
% Written for single time-point images. 

% 7-23-18: changed the filtering step. Using diffusion in 3D code. Also,
% calculating gaussian gradient in 3D using separable filters. Now we don't
% have to loop over z. Hopefully this is a bit faster! ALSO, changing the
% thresholding step to use an absolute value rather than fraction of total
% luminance. 
% Also, keeping reader open during loop over time. Might be taking a long
% time to open and close image reader. 

function obj = segment_3D_gpu(obj, params, step)



%% Pre-processing using bio-formats. 

%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = bfGetReader(obj.exp_info.img_file);
series = 1;


%Get the image size of this series. 

T = reader.getSizeT;

 %Get step.debug
if(step.debug)
    params
    %T = 1;
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
            inner_function( exp_info, Z_cell{t}, params, step, t, reader)
        end
    else
        disp('New frames to calculate:')
        disp(num2str(new_ts))
        for t = new_ts
            inner_function( exp_info, Z_cell{t},params,step,t, reader)
        end
    end

%Close reader at very end. 
reader.close();
end

%% Inner function to run 
function inner_function( exp_info, planes, params, step, t, reader)
    
%Display 
disp(['Started frame: ',num2str(t)])

%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%  
         
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
        %Mean filtering 
        MeanFilterSensitivity = params.MeanFilterSensitivity;
        MeanFilterNeighborhood = params.MeanFilterNeighborhood;

        %Imclose
        imclose_r = params.imclose_r;
        
        xscale         = exp_info.pixel_size(1);
        yscale         = exp_info.pixel_size(2);
        zscale         = exp_info.pixel_size(3);

%%%%%%%%%%%%%%% START PROCESSING %%%%%%%%%%%%%%%%
    image_bits     = reader.getBitsPerPixel;
    size_x = reader.getSizeX;
    size_y = reader.getSizeY;
    ZSlicesinStack = length(planes);

    %Get image planes for this time point. 
    I = zeros(size_y,size_x,ZSlicesinStack);

    %Get the bio-formats image index corresponding to this z-stack:
    for i = 1:ZSlicesinStack
        this_plane_img = bfGetPlane(reader,planes(i));
        %If there are zeros in this plane due to weird stitching, then
        %adjust. 
        sel = this_plane_img == 0;
        if(sum(sel(:))>0)
            this_plane_img(sel) = mean( this_plane_img(~sel) );
        end
        
        I(:,:,i)     = this_plane_img;
    end
    
    %Reduce i for debugging
    if(step.debug)
       % I = I(1:600,1:600,:);
    end
    Ie = zeros(size(I));
     
    
    %% PRE-PROCESSING STEPS. Filters, and removing outliers. 
       
 
    
    %PRE-PRfiOCESS image with mean filter if desired
    if(step.MeanFilter)
        tmp = adaptthresh(uint16(I),MeanFilterSensitivity,'NeighborhoodSize',MeanFilterNeighborhood);
        
        %%%% This step might help bring everything to 16bit levels. Try
        %%%% manually changing bit size to 16 here. 
        I = tmp.*65000;
        image_bits = 16;
    end

    %CLAHE
    if(step.CLAHE)
        I_sm = zeros(size(I));
        for i = 1:ZSlicesinStack
            this_plane = I(:,:,i);
            norm_max = max(this_plane(:));
            norm_min = min(this_plane(:));
            
            norm_plane = ( this_plane - norm_min ) ./ (norm_max - norm_min );
            tiles = round([size_y / params.tile_size, size_x / params.tile_size]);
            I_sm(:,:,i) = adapthisteq(norm_plane,'NumTiles',tiles,'ClipLimit',params.ClipLimit); %,'Distribution','exponential','Alpha',0.1);
            
        end
        
        I = I_sm.*65535;
        
        
    end
    
    if(step.subtract_bg)
        sel = I < params.bg;
        I(sel) = 0;
    end
    
    disp('start filtering')
    
    %Start timing.
    tic
    
    %Speed up using gpuArrays. 
    I_gpu = gpuArray(I);
    G = gpuArray(Ie);
    
    %Median filtering and diffusion in 2D. 
    for i  = 1:ZSlicesinStack
        G(:,:,i) = medfilt2(I_gpu(:,:,i),params.median_filter(1:2));
    end
    G = gather(G);
    Diff_im = gpuArray(Ie);
    for i =  1:ZSlicesinStack
        Diff_im(:,:,i) = diffusioncode(I_gpu(:,:,i), diffuse_iterations, 0.1429, kappa1, kappa2, option);
    end
    %Clear out for memory
    clear I_gpu
    Diff_im = gather(Diff_im);
    
    %Rescale to [0,1]
    Fim = mat2gray(Diff_im);
    
    %Get gauss gradient. Using separable filters. sigmagradient now 3 component vector
    [imx,imy]=gaussgradient2D_sep(Fim,sigmagradient);
    
    %Total magnitude. 
    Mag_fim = hypot(imx,imy);
    
    %Laplacian. 
    [L_imxx, L_imxy] = gaussgradient2D_sep(imx,sigmagradient);
    [L_imyx, L_imyy] = gaussgradient2D_sep(imy,sigmagradient);

    
    Mag_Lim = L_imxx + L_imyy;% + L_imxy + L_imyx; % + L_imzz;
    Mag_Lim = Mag_Lim.*(Mag_Lim>0);
    
    Det_hessian = L_imxx.*L_imyy - L_imxy.*L_imyx;
    Det_hessian = Det_hessian.*(Det_hessian<0);

    X=gamma-((alpha*Mag_fim + beta*Mag_Lim + epsilon*abs(Det_hessian)) ./ delta);
    Multi=0.5*(tanh(X)+1);% masking function
    
    %Apply mask and gather back to cpu.
    J=double(G).*Multi; % masked image applied on smoothened image
    
    %Simple binarization. 
    BW  = J >= thrshlevel;
 
    %Collect stats. 
    stats = regionprops(BW,'Centroid','Area','PixelList','PixelIdxList');
    
    %Remove volumes below size threshold. 
    volumes = [stats.Area]';
    sel = volumes >= AbsMinVol;
    stats = stats(sel);
    
    %Rebuild BW
    BW = zeros(size(BW));
    idx = cat(1,stats.PixelIdxList);
    BW(idx) = 1;
    BW = logical(BW);
    
    if(step.imclose)
        %Now use imclose to remove gaps. Operate on individual cells so
        %that we minimize connecting cells that are close together. 
        se= strel('disk',imclose_r);
        new_BW= gpuArray(BW);
        
        stats = regionprops(gather(new_BW),'Centroid','Area','PixelList','PixelIdxList');
        
        for c = 1:length(stats)
            X = stats(c).PixelList(:,1);
            Y = stats(c).PixelList(:,2);
            Z = stats(c).PixelList(:,3);

            x_range = [min(X):max(X)];
            y_range = [min(Y):max(Y)];
            z_range = [min(Z):max(Z)];

            %Perform morphological closing.
            new_BW(y_range,x_range,z_range) = imclose(new_BW(y_range,x_range,z_range),se);
        end
        
    else
        new_BW = gpuArray(BW);
            
    end

    %Now fill in some holes. 
    for i = 1:1:ZSlicesinStack
        %Also fill holes laterally. Vertical holes can be
        %problematic if we just to imfill with 3D stack
        new_BW(:,:,i) = imfill(new_BW(:,:,i),'holes');
    end
    
    BW = gather(new_BW);

    %End timing
    toc
    
    %Clear out GPU
    gpuDevice();
    clear I_gpu
  
        
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
    
    %New volume list
    volumes = [stats.Area]';   
    liar_list = find( volumes > AbsMaxVol);
    
    %% Apply 3D watershedding. 
    if( step.watershed3D)
        disp('Starting watershed')
        %Make bw of clusters only
        bw_clust = zeros(size(BW));
        idx = cat(1,stats(liar_list).PixelIdxList);
        bw_clust(idx) = 1;
        %Anisotropic distance transform
        scales = [exp_info.pixel_size(1:2),exp_info.pixel_size(3)*params.z_effect];
        imgDist = -bwdistsc_gpu(~bw_clust,scales);

        %Smoothing. Med filter is actually important here! 
        imgDist = medfilt3(imgDist,[5,5,5]);
        
        %Get seeds  %%ORIGINAL params were 0.7,6. With medfilt3 = 5,5,5. 
        mask = imextendedmin(imgDist,params.h_min_depth,params.h_min_conn); %Seems like smaller neighborhood works better?
        imgDist = imimposemin(imgDist,mask);
        imgDist(~bw_clust) = -inf;
        %Perform marked watershed
        imgLabel = watershed(imgDist);
        %Replace background again...
        imgLabel(~bw_clust) = NaN;
        stats_fused = regionprops(imgLabel,'Centroid','Area','PixelList','PixelIdxList');
        stats_fused = stats_fused(2:end);
        
        %Now get rid of original regions that were fused. 
        keep = setdiff([1:length(stats)],liar_list);
        stats = stats(keep);
        %Now append the fused nuclei stats. 
        stats = [stats;stats_fused];
        %Now we need to adjust BW... 
            BW = zeros(size(BW));
            %Add in new regions
            new_idx = cat(1,stats.PixelIdxList);
            BW(new_idx) = 1;
        disp('Finished watershed')
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
    rg_all = [];
    %Empty BW. 
    BW = zeros(size(BW));
    for i = 1:length(stats)
            if (length(stats(i).PixelIdxList) < AbsMinVol)
                continue
            end

            %Check radius of gyration.
            [y,x,z] = ind2sub(size(BW),stats(i).PixelIdxList);
            %Calculate radius of gryation
            pos=[x,y,z];
            com=[mean(x),mean(y),mean(z)];
            rg = sqrt(1/length(y).*sum(sum( (pos-com).^2,2)));
            rg_all = [rg_all, rg];
            if(step.rg_threshold & rg > params.rg_threshold)
                continue
            end
            
            %Fill in BW. 
            BW(stats(i).PixelIdxList) = 1;
            frame_obj.PixelIdxList{counter} = stats(i).PixelIdxList;
            frame_obj.centroids{counter}    = stats(i).Centroid;
            counter = counter + 1;
    end

        %Add final binarized image to frame_obj for safe keeping.
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

try
    save(fname,'frame_obj');
catch
    %For large files...
    save(fname, 'frame_obj','-v7.3');
end
end

