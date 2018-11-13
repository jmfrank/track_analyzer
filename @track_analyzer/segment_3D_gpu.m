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

% 9-3-18: Major speed up by performing imclose and watershedding blobs individually. Also, fixed a
% bug with the imclose section: was copying the region of the blob instead
% of the individual voxels which sometimes included regions of other
% blobs...
%
% 9-3-18: We should write this function to just take an image stack and use
% a different wrapper. This would im prove usability. 

function obj = segment_3D_gpu(obj, params, step, FORCE_FRAMES)

params = default_params(params);
step   = default_step(step);

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
    ZSlicesinStack = reader.getSizeT;

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
        
    
%% Loop over images. 

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

%If user specifies list of frames as last input, then for frames accordinly
if nargin==4
    new_ts = FORCE_FRAMES;
end


%Now run loops
disp('New frames to calculate:')
disp(num2str(new_ts))
for t = new_ts
    inner_function( exp_info, Z_cell{t},params,step,t, reader)
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
        plane_mean(i) = mean(this_plane_img(:));
        plane_std(i)  = std(double(this_plane_img(:)));
        I(:,:,i)     = this_plane_img;
    end
    
     
    
    %% PRE-PROCESSING STEPS. Filters, and removing outliers. 

    %PRE-PROCESS image with mean filter if desired
    if(step.MeanFilter)
        
        if image_bits==8
            if ~strcmp(class(I),'uint8')
                I = uint8(I);
            end
            tmp = adaptthresh(I,MeanFilterSensitivity,'NeighborhoodSize',MeanFilterNeighborhood);
        elseif image_bits==16
            tmp = adaptthresh(uint16(I),MeanFilterSensitivity,'NeighborhoodSize',MeanFilterNeighborhood);
        end
        
        %%%% This step might help bring everything to 16bit levels. Try
        %%%% manually changing bit size to 16 here. 
        I = tmp.*65000;
        image_bits = 16;
    end

    %CLAHE
    if(step.CLAHE)
        I_sm = zeros(size(I));
        %Weights. 
        lower_bound = min(plane_mean)-0.5.*min(plane_std);
        max_bound   = max(plane_mean);
        weights = (plane_mean-lower_bound) ./ (max_bound-lower_bound);
        for i = 1:ZSlicesinStack
            this_plane = I(:,:,i);
            norm_max = max(this_plane(:));
            norm_min = min(this_plane(:));
            
            norm_plane = ( this_plane - norm_min ) ./ (norm_max - norm_min );
            tiles = round([size_y / params.tile_size, size_x / params.tile_size]);
            I_sm(:,:,i) = adapthisteq(norm_plane,'NumTiles',tiles,'ClipLimit',params.ClipLimit).*weights(i); %,'Distribution','exponential','Alpha',0.1);
            
        end
        
        I = I_sm.*65535;
       
    end
    
    if(step.subtract_bg)
        sel = I < params.bg;
        I(sel) = 0;
    end
    
    disp('start filtering')
    tic
    
    %Speed up using gpuArrays. 
    I_gpu = gpuArray(I);
    Ie = zeros(size(I));
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
    
    %Estimate threshold by histogram. 
    if step.threshold_by_histogram 
        P = prctile(J(:),params.percentile);
        thrshlevel = P;
    end
    
    %Simple binarization. 
    BW  = J >= thrshlevel;
 
    if step.debug
        
       imshow3D(J) 
    end
    
    
    %Collect stats. 
    stats = regionprops(BW,'Centroid','Area','PixelList','PixelIdxList');
    
    %Remove VERY small volumes
    volumes = [stats.Area]';
    sel = volumes >= AbsMinVol/4;
    stats = stats(sel);
    
    %Rebuild BW
    BW = zeros(size(BW));
    idx = cat(1,stats.PixelIdxList);
    BW(idx) = 1;
    BW = logical(BW);
    
    if(step.imclose)
        %Now use imclose to remove gaps. Operate on individual cells so
        %that we minimize connecting cells that are close together. 

        %Ellipsoidal structuring element,. 
        %se = strel3d_ellipsoid(imclose_r,imclose_r,round(imclose_r/3));
        se = strel('disk',imclose_r,8);
        
        new_BW= zeros(size(BW));
        %Re-dundant regionprops? 
        stats = regionprops(BW,'Centroid','Area','PixelList','PixelIdxList');

        for c = 1:length(stats)
            X = stats(c).PixelList(:,1);
            Y = stats(c).PixelList(:,2);
            Z = stats(c).PixelList(:,3);

            x_range = [min(X):max(X)];
            y_range = [min(Y):max(Y)];
            z_range = [min(Z):max(Z)];
            
            %Need to make a sub-img! 
            sub_img = false(length(y_range),length(x_range),length(z_range));
            %Shift og pixels to sub_img pixels. 
            X = X - x_range(1) + 1;
            Y = Y - y_range(1) + 1;
            Z = Z - z_range(1) + 1;
            %Get index. 
            ind = sub2ind(size(sub_img),Y,X,Z);
            sub_img(ind) = 1;
            %Close. 
            sub_img = imclose(sub_img,se);
            
            %Find closed coordinates. 
            [Y,X,Z] = ind2sub(size(sub_img),find(sub_img));
            %Shift back into place. 
            X = X + x_range(1) - 1;
            Y = Y + y_range(1) - 1;
            Z = Z + z_range(1) - 1;
            %New index. 
            ind = sub2ind(size(BW),Y,X,Z);
            %Fill in the large image. 
            new_BW(ind) = 1;
        end
        new_BW = gpuArray(new_BW);
        
    else
        new_BW = gpuArray(BW);
            
    end

    %Now fill in some holes. 
    for i = 1:1:ZSlicesinStack
        %Also fill holes laterally. Vertical holes can be
        %problematic if we just to imfill with 3D stack
        new_BW(:,:,i) = imfill(new_BW(:,:,i),'holes');
    end
    
    %Collect stats again. 
    stats = regionprops(logical(gather(new_BW)),'Centroid','Area','PixelList','PixelIdxList');
    %Remove VERY small volumes
    volumes = [stats.Area]';
    sel = volumes >= AbsMinVol;
    stats = stats(sel);
    
    %Rebuild BW
    BW = false(size(BW));
    BW( cat(1,stats.PixelIdxList) ) = 1;

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
    
    %Plot liar tex
    if(step.debug)
        plot_liar_text( stats, liar_list );
    end
    
    %% Apply 3D watershedding. Much faster doing individually for each object <8-31-18. JMF. 
    if( step.watershed)
        disp('Starting watershed')
        tic
        
        %Distance transform scales. 
        scales = [exp_info.pixel_size(1),exp_info.pixel_size(2),exp_info.pixel_size(3)*params.z_effect];

        %Counter of new objects found by watershedding
        new_objects = 0;

        %Empty fused object imglabel
        fused_imgLabel = uint8(zeros(size(BW))); 

        %Loop over liars
        for i = 1:length(liar_list)

            %This liar idx
            idx = liar_list(i);

            %Sub_img containing this liar_volume. 
            Px = stats( idx ).PixelList;

            %Pixel ranges. 
            min_range = min(Px,[],1);
            max_range = max(Px,[],1);

            x_range = min_range(1):max_range(1);
            y_range = min_range(2):max_range(2);
            z_range = min_range(3):max_range(3);

            %Create a sub image of BW. 
            sub_bw = zeros(length(y_range),length(x_range),length(z_range));
            
            %Add in the ones from this blob. 
            X = Px(:,1) - x_range(1) + 1;
            Y = Px(:,2) - y_range(1) + 1;
            Z = Px(:,3) - z_range(1) + 1;
            ind = sub2ind(size(sub_bw),Y,X,Z);
            sub_bw(ind) = 1;
            
            %Distance transform (faster without using gpu when image is
            %small. 
            imgDist = -bwdistsc(~sub_bw,scales);

            %Smoothing. Med filter is actually important here! 
            imgDist = medfilt3(imgDist,[3,3,3]);    

            %Get seeds  %%ORIGINAL params were 0.7,6. With medfilt3 = 5,5,5. 
            mask = imextendedmin(imgDist,params.h_min_depth,params.h_min_conn); %Seems like smaller neighborhood works better?
            imgDist = imimposemin(imgDist,mask);
            imgDist(~sub_bw) = -inf;

            
            %Perform marked watershed
            imgLabel = watershed(imgDist);
            %Replace background again?
            imgLabel(~sub_bw) = NaN;

            %Putting objects into new compiled BW. Figure out how many
            %objects exist after watershedding. Infer from the values of
            %imgLabel. imglabel=0 or 1 is background or envelope. 

            %Number of objects
            num_objects = double(max(imgLabel(:))-1 ); %Subtract 1 cause 0 and 1 are nothing. 

            %Ignore non-object pixels. 
            sel = imgLabel > 1;
            imgLabel(~sel) = 0;

            %Shift values by the on going counter. 
            imgLabel(sel) = imgLabel(sel) + new_objects;

            %Now place the imgLabel into the fused BW img back where it
            %came from. 
            fused_imgLabel(y_range, x_range, z_range ) = fused_imgLabel(y_range, x_range, z_range ) + imgLabel;

            %Adjust counter. 
            new_objects = new_objects + num_objects;
        end
        
        %Keep non liar stats
        keep = setdiff([1:length(stats)],liar_list);
        stats = stats(keep); 

        %Collect stats on watershedded regions. 
        stats_fused = regionprops(fused_imgLabel,'Centroid','Area','PixelList','PixelIdxList');
        stats_fused = stats_fused(2:end); %Ignore first entry. 
        %Append to stats. 
        stats = [stats; stats_fused];
        %Now we need to adjust BW... 
        BW = zeros(size(BW));
        %Add in new regions
        new_idx = cat(1,stats.PixelIdxList);
        BW(new_idx) = 1;
        toc        
        
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
        if exist('disp_str','var'); clearString(disp_str); end
        disp_str = ['Finished frame:     ',num2str(t)];
        disp(disp_str)
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


%% Plotting text over image
function    plot_liar_text( stats, liar_list )

delete(findobj(gca,'Type','text'))
hold on
for i = 1:length(liar_list)
    ctr = stats( liar_list(i) ).Centroid;
    text(ctr(1),ctr(2),num2str( liar_list(i)),'color','g')
end
hold off

end