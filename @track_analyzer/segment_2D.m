% 2D Segmentation of cell nuclei. Based on: 'Object segmentation and ground
% truth in 3D embryonic imaging'. 

% This version is written for 2D images. If there's a z-stack, then just
% max-int project first. 

% 10-26-18: now takes vector of threshold values for each time frame.
function obj = segment_2D(obj, params, step, FORCE_FRAMES)

params=default_params(params);
step=default_step(step);

%% Pre-processing using bio-formats. 

%Generate reader. Use memoizer. Oif files can be really slow. Since segmentation is
%the first part of processing, there should be memo for all subsequent
%steps. 
bfInitLogging;

if step.use_specific_img    
    loci.common.DebugTools.setRootLevel('WARN');
    reader = bfGetReader();
    reader = loci.formats.Memoizer(reader,0);
    reader.setId( params.img_file );
    series = 1;    
else
    loci.common.DebugTools.setRootLevel('WARN');
    reader = bfGetReader();
    reader = loci.formats.Memoizer(reader,0);
    reader.setId( obj.exp_info.img_file );
    series = 1;
end




%Get the image size of this series. 

T = reader.getSizeT;

%% Generate im_info structure
    ZSlicesinStack = reader.getSizeZ;
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

%If user specifies list of frames as last input, then for frames accordinly
if nargin==4
    new_ts = FORCE_FRAMES;
    
    if any( new_ts > T )
        error(['frame requested exceeds image frame count: ',num2str(T)]);
    end
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
            inner_function( exp_info, Z_cell{t},params,step,t,reader)
        end
    end
disp(['Max frames: ',num2str(T)]);

%Add params to exp_info.
obj.exp_info.seg_params=params;
obj.exp_info.seg_steps=step;
%Always clear out flags in after new segmentation. 
obj = obj.clear_flags;
%Auto-save. 
obj.save; 

obj.save;


end

%% Inner function to run 
function inner_function( exp_info, planes, params, step, t,reader)
    
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

        %%Perona-Malik based non-linear isotropic diffusion filter
        diffuse_iterations=params.diffuse_iterations;
        kappa1=params.kappa1;
        kappa2=params.kappa2;
        option=params.option;

        % nconn value for neighborhood calculations
        nconn_BW = params.nconn_BW;

        MeanFilterSensitivity = params.MeanFilterSensitivity;
        MeanFilterNeighborhood = params.MeanFilterNeighborhood;

        imclose_r = params.imclose_r;   

%%%%%%%%%%%%%%% START PROCESSING %%%%%%%%%%%%%%%%
    image_bits     = reader.getBitsPerPixel;
    size_x = reader.getSizeX;
    size_y = reader.getSizeY;
    ZSlicesinStack = length(planes);

    %Get image planes for this time point. 
    I = zeros(size_y,size_x,ZSlicesinStack);

    %Get the bio-formats image index corresponding to this z-stack:
    if ZSlicesinStack > 1
        
        for i = 1:ZSlicesinStack
            this_plane_img = bfGetPlane(reader,planes(i));
            I(:,:,i)     = this_plane_img;
        end
        %Project z-dimension. 
        I = max(I,[],3);
    else
        I = bfGetPlane(reader,planes);
    end
    
    
    %Reduce i for debugging
    if(step.debug)
        %I = I(1:300,1:300,:);
    end
    

    
    Ie = zeros(size(I));
    G = Ie; Diff_im = Ie; Fim = Ie; Mag_fim = Ie; Mag_Lim = Ie; Det_hessian = Ie;
    Sum = Ie; Multi = Ie; J = Ie;
     
    
    %% PRE-PROCESSING STEPS. Filters, and removing outliers. 
       
    %PRE-PROCESS image with mean filter if desired
    if(step.MeanFilter)
        
        if image_bits==8
            if ~strcmp(class(I),'uint8')
                I = uint8(I);
            end
            tmp = adaptthresh(I,MeanFilterSensitivity,'NeighborhoodSize',MeanFilterNeighborhood(1:2));
        elseif image_bits==12
            
            %Rescale to 16 bit? 
            I = uint16(I* (2^16-1)/(2^12-1));
            
        elseif image_bits==16
            if ~strcmp(class(I),'uint16')
                I = uint16(I);
            end
        end
        
        tmp = adaptthresh(I,MeanFilterSensitivity,'NeighborhoodSize',MeanFilterNeighborhood(1:2));

        %%%% This step might help bring everything to 16bit levels. Try
        %%%% manually changing bit size to 16 here. 
        I = tmp.*65000;
    end
    
    %Subtract a background value. 
    if(step.subtract_bg)
       sel = I <= params.bg;
       I(sel) = 0;
    end
    
    %CLAHE
    if(step.CLAHE)
        I_sm = zeros(size(I));
        mu = mean(I(:));
        i_std = std(I(:));
        norm_min = min(I(:));
        norm_max = max(I(:));
        %Set the max val to 4std above mean?
        %norm_max = 4*i_std + mu;
            
        norm_plane = ( I - norm_min ) ./ (norm_max - norm_min );
        tiles = round([size_y / params.tile_size, size_x / params.tile_size]);
        I_sm = adapthisteq(norm_plane,'NumTiles',tiles,'ClipLimit',params.ClipLimit); %,'Distribution','exponential','Alpha',0.1);
            
        %Replace background 
        if(step.subtract_bg)
            I_sm(sel) = 0;
        end
       
        I = I_sm.*65535;
       
    end   
    
    %% Diffusion 

    %%Step1: Denoising Filters
    %G(:,:,zslice) = imfilter(I(:,:,zslice), h,'replicate');
    %%Optional filters median or deconvlucy
    G = medfilt2(I, params.median_filter(1:2));
    %G(:,:,zslice) = deconvlucy(G(:,:,zslice),h); %%use the same 'h' as for Gaussian or define new variable

    %%Step2: Perona & Malik non-linear isotropic diffusion filter, refer to
    %%diffusioncode.m for details, here alter only diffuse_iterations & kappa1       
    Diff_im = diffusioncode(I, diffuse_iterations, 0.1429, kappa1, kappa2, option);

    %%Step3a: Gauss gradient, 1st derivative of the image, refer to gaussgradient.m for details
    Fim=mat2gray(Diff_im);
    [imx,imy]=gaussgradient(Fim,sigmagradient(1));
    Mag_fim = hypot(imx,imy);

    %%Step3b: Laplacian
    [L_imxx, L_imxy] = gaussgradient(imx,sigmagradient(1));
    [L_imyx, L_imyy] = gaussgradient(imy,sigmagradient(1));
    Mag_Lim = L_imxx+L_imyy;  %Laplacian (trace of the 2D matrix)
    Mag_Lim=Mag_Lim.*(Mag_Lim>0);

    %%Step3c: Hessian Determniant
    Det_hessian=L_imxx.*L_imyy-L_imxy.*L_imyx;
    Det_hessian=Det_hessian.*(Det_hessian<0);
    %%Step4: Masking function using tanh on the Summed Derivatives from Step3
    X=gamma-((alpha*Mag_fim + beta*Mag_Lim + epsilon*abs(Det_hessian)) /delta);
    Multi=0.5*(tanh(X)+1);% masking function
    J=(double(G).*Multi); % masked image applied on smoothened image

    if(step.debug)
        figure(10)
        if isfield(params,'percentile')
            try
                K=params.percentile(t);
            catch
                K=params.percentile;
            end
        else
            K=90;
        end
        imshow3D_filter(J,prctile(J(:),K));
        
    end 

    %% Different ways to determine threshold. 
    
    %Iterative local thresholding. 
    if step.iterative_thresholding
        
        if length(params.percentile)>1
            p_val = params.percentile(t);
        else
            p_val = params.percentile;
        end

        %Enforce more criteria for thresholing cells. 
        real_bg = mean(J(sel));
        %Threshold starting point. 
        params.thresh_start = prctile(J(:),p_val);
        %Incremental change in threshold. 
        params.increment = params.thresh_start*0.01;
        %First smooth image. 
        I_sm = imgaussfilt(I,params.I_sm_sigma);
        %Perform iterative thresholding. 
        BW = iterative_thresholding(I_sm, J,real_bg, params );
        
        
    %Histogram thresholding. 
    elseif step.threshold_by_histogram 
        
       
        	if length(params.percentile)>1
                p_val = params.percentile(t);
            else
                p_val = params.percentile;
            end
            P = prctile(J(:),p_val);
            thrshlevel = P;
            BW  = J >= thrshlevel;
    
    %Simple binarization. 
    else
        if length(params.thrshlevel)>1
            thrshlevel=params.thrshlevel(t);
        else
            thrshlevel=params.thrshlevel;
        end
        BW  = J >= thrshlevel;
    end
    
    %% smoothing, filling holes, removing bad blobs. 
          
    %Collect stats. 
    stats = regionprops(BW,'Centroid','Area','PixelList','PixelIdxList');

    %Remove very low volumes. 
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
        se= strel('disk',imclose_r,8);
        new_BW= zeros(size(BW));


        for c = 1:length(stats)
            X = stats(c).PixelList(:,1);
            Y = stats(c).PixelList(:,2);

            x_range = [min(X):max(X)];
            y_range = [min(Y):max(Y)];

            %Need to make a sub-img! 
            sub_img = false(length(y_range),length(x_range));
            %Shift og pixels to sub_img pixels. 
            X = X - x_range(1) + 1;
            Y = Y - y_range(1) + 1;   
            %Get index. 
            ind = sub2ind(size(sub_img),Y,X);
            sub_img(ind) = 1;
            %Close. 
            sub_img = imclose(sub_img,se);
            
            %Find closed coordinates. 
            [Y,X] = ind2sub(size(sub_img),find(sub_img));
            %Shift back into place. 
            X = X + x_range(1) - 1;
            Y = Y + y_range(1) - 1;
            
            %New index. 
            ind = sub2ind(size(BW),Y,X);
            
            %Fill in the large image. 
            new_BW(ind) = 1;            
            
        end
        BW = logical(new_BW);
    end

    %Now fill in the holes. 
    BW = imfill(BW,'holes');

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
    
    %% Looking for high intensity outliers. 
    if(step.FilterOutliers)
       
        %Estimate background value. 
        bg_img = I(~BW);
        mean_bg_value = mean(bg_img(:));
        
        %Find outliers using a high threshold. 
        high = I > params.OutlierThreshold;
        
        %Dilate 
        se = strel('rectangle',params.MeanFilterNeighborhood );
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
    
    %% Connect regions, find objects that are are big enough to be doubles. 

    Con=bwconncomp(logical(BW),nconn_BW);
    stats = regionprops(Con,'Centroid','Area','PixelList','PixelIdxList');
    numobj=numel(stats); %total number of segmented objects
 
    %New volume list
    volumes = [stats.Area]';   
    liar_list = find( volumes > AbsMaxVol);
    
   %% Apply 2D watershedding. This only looks at the blobs that are above the AbsMaxVol threshold. 
    if( step.watershed)
        disp('Starting watershed')
        %Make bw of clusters only
        bw_clust = zeros(size(BW));
        idx = cat(1,stats(liar_list).PixelIdxList);
        bw_clust(idx) = 1;
        %distance transform
        imgDist = -bwdist(~bw_clust);
        %Smoothing
        imgDist = medfilt2(imgDist,[3,3]);
        %Get seeds 
        mask = imextendedmin(imgDist,params.h_min_depth,params.h_min_conn); %Seems like smaller neighborhood works better?
        imgDist = imimposemin(imgDist,mask);
        imgDist(~bw_clust) = -inf;
        %Perform marked watershed
        imgLabel = watershed(imgDist);
        %Remove background again.
        imgLabel(~bw_clust) = NaN;
        stats_fused = regionprops(imgLabel,'Centroid','Area','PixelList','PixelIdxList');
        %skip first entry. 
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
    
    %% Error if there's no cells.
    if(~isfield(stats,'PixelList'))
        error('No cells found!')
    end

    %% Merge blobs that are very close together. 
    if step.merger
        
        %Centroids. 
        all_ctrs = cat(1,stats.Centroid); 
        
        %Pairwise distance. 
        D = pdist2(all_ctrs,all_ctrs);
        
        %Select pair-wise distances less than threshold. 
        S = D < params.merge_dist;
        
        %Get pairs. 
        [a,b] = find( tril( S ) );
        
        %Entries to get rid of. 
        bad_list = false(1,length(stats));
        c=0;
        
        %Merge algorithm. Loop over pairs, merge if result isn't too large.
        for m = 1:length(a)
            id_a = a(m);
            id_b = b(m);
            w_a = length(stats(id_a).PixelIdxList);
            w_b = length(stats(id_b).PixelIdxList);
            merge_vol = w_a + w_b;
            
            if merg_vol < AbsMaxVol
                bad_list(id_a) = 1;
                bad_list(id_b) = 1;
                new_list = [stats(id_a).PixelIdxList; stats(id_b).PixelIdxList];
                new_centroid = w_a/merge_vol*stats(id_a).Centroid + w_b/merge_vol*stats(id_b).Centroid;
                new_pixel_list = [stats(id_a).PixelList; stats(id_b).PixelList];
                c=c+1;
                new_stats(c).Centroid=new_centroid;
                new_stats(c).PixelList = new_pixel_list;
                new_stats(c).PixelIdxList = new_list;
            end
        end
        
        %Now remove bad objects, and append new ones. 
        stats = stats(~bad_list);
        stats = [stats; new_stats];
        
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
            [y,x] = ind2sub(size(BW),stats(i).PixelIdxList);
            %Calculate radius of gryation
            pos=[x,y];
            com=[mean(x),mean(y)];
            
            rg = sqrt(1/length(y).*sum(sum( (pos-com).^2,2)));
            rg_all = [rg_all,rg];
            
            if(step.rg_threshold && rg > params.rg_threshold)
                continue
            end
            
            %Fill in BW. 
            BW(stats(i).PixelIdxList) = 1;
            frame_obj.PixelIdxList{counter} = stats(i).PixelIdxList;
            frame_obj.centroids{counter}    = stats(i).Centroid;
            
            counter = counter + 1;
    end
        %Trace boundaries. 
        C = bwboundaries(BW);
        if ~isempty(C)
            C = cellfun(@(x) [smooth(x(:,2)),smooth(x(:,1))],C,'uniformoutput',0);
            c_ctr = cellfun(@(x) [mean(x(:,1)),mean(x(:,2))],C,'uniformoutput',0);
            %Match centroids. 
            D = pdist2(cat(1,c_ctr{:}),cat(1,frame_obj.centroids{:}));
            [~,idx] = min(D);
            frame_obj.contours=C( idx );
        end
        
        
        
        %Add final binarized image to frame_obj for save keeping
        frame_obj.BW = BW;

        %Save frame_obj as done before. 
        fname = ['frame_',sprintf('%04d',t),'.mat'];
        if(~exist(exp_info.nuc_seg_dir,'dir'))
            mkdir(exp_info.nuc_seg_dir)
        end
        parsave([exp_info.nuc_seg_dir,fname],frame_obj)
        if isempty(stats)
            nFound=0;
        else
            nFound=length(frame_obj.contours);
        end
        disp(['Finished frame:     ',num2str(t),' Found ',num2str(nFound),' cells.'])

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


%Plotting text over image
function plot_liar_text( stats, liar_list )

delete(findobj(gca,'Type','text'))
hold on
for i = 1:length(liar_list)
    ctr = stats( liar_list(i) ).Centroid;
    text(ctr(1),ctr(2),num2str( liar_list(i)),'color','g')
end
hold off

end
