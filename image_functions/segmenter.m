% class for performing segmentation steps. 
classdef segmenter < handle  
    
    properties (SetAccess = private)
        
        %BW
        
        %img
        
        %filtered
        
        stats
        
        image_bits
        
        t % frame.
        
        mask
        
        steps %< internal tracking of all steps called
    end
    
    properties (SetAccess = public)
       
        %changeable parameters. 
        params 
        
        msk_pxs

        msk_centroids
        
        BW
        
        img
        
        filtered  

        seg_type

        spot_fits
        
    end
    
    
    
    methods (Access=public)
        
        % Constructer
        function obj = segmenter( I, image_bits, params, t, seg_type, msk_pxs  ) 
            
            obj.img = I;
            obj.image_bits= image_bits;
            obj.params=params;
            obj.t=t;
            obj.seg_type=seg_type;
            obj.BW = false(size(obj.img));
            %Default mask is True! Hope this doens't screw things up...  
            obj.mask = true(size(obj.img));    
            if nargin > 5 % need to provide msk_pxs. 
                obj.add_mask(msk_pxs);  
            end
        end
        
        % Image processing functions. 
        function mean_filter(obj)
            
            % Check bit-depth. 
            max_bits = max(obj.img(:));
            
            if obj.image_bits==8
                if ~strcmp(class(obj.img),'uint8')
                    obj.img = uint8(obj.img);
                end
                tmp = adaptthresh(obj.img,obj.params.MeanFilterSensitivity,'NeighborhoodSize',obj.params.MeanFilterNeighborhood);
            elseif obj.image_bits==12

                %Rescale to 16 bit because matlab doesn't have 12-bit
                %object. 
                obj.img = uint16(obj.img* (2^16-1)/(2^12-1));
                tmp = adaptthresh(obj.img,obj.params.MeanFilterSensitivity,'NeighborhoodSize',obj.params.MeanFilterNeighborhood);

            elseif obj.image_bits==16
                if ~strcmp(class(obj.img),'uint16')
                    obj.img = uint16(obj.img);
                end
                tmp = adaptthresh(obj.img,obj.params.MeanFilterSensitivity,'NeighborhoodSize',obj.params.MeanFilterNeighborhood);

            end

            %%%% This step might help bring everything to 16bit levels. Try
            %%%% manually changing bit size to 16 here. 
            obj.img = tmp.*65000;            
            
        end
        
        function subtract_background(obj)
        
            sel = obj.img <= obj.params.bg;
            obj.img(sel) = 0;          
            
        end
        
        function clahe(obj)
            
             % error if 3D. 
            if length(size(obj.img))==3

                error('CLAHE function not support in 3D');

            end

            I_sm = zeros(size(obj.img));
            mu = mean(obj.img(:));
            i_std = std(obj.img(:));
            norm_min = min(obj.img(:));
            norm_max = max(obj.img(:));
            %Set the max val to 4std above mean?
            %norm_max = 4*i_std + mu;
            [size_y,size_x,size_z] = size(obj.img);
            norm_plane = ( obj.img - norm_min ) ./ (norm_max - norm_min );
            tiles = round([size_y /obj.params.tile_size, size_x / obj.params.tile_size]);
            I_sm = adapthisteq(norm_plane,'NumTiles',tiles,'ClipLimit',obj.params.ClipLimit); %,'Distribution','exponential','Alpha',0.1);

            obj.img = I_sm.*65535;
            
        end
        
        function diff_gradients(obj)
            
            % Check dimensionality. 
            dims = length(size(obj.img));

            switch dims

                case 2

                    G = medfilt2(obj.img, obj.params.median_filter(1:2));
                    %G(:,:,zslice) = deconvlucy(G(:,:,zslice),h); %%use the same 'h' as for Gaussian or define new variable

                    %%Step2: Perona & Malik non-linear isotropic diffusion filter, refer to
                    %%diffusioncode.m for details, here alter only diffuse_iterations & kappa1       
                    Diff_im = diffusioncode(obj.img, obj.params.diffuse_iterations, 0.1429, obj.params.kappa1, obj.params.kappa2, obj.params.option);

                    %%Step3a: Gauss gradient, 1st derivative of the image, refer to gaussgradient.m for details
                    Fim=mat2gray(Diff_im);
                    [imx,imy]=gaussgradient(Fim,obj.params.sigmagradient(1));
                    Mag_fim = hypot(imx,imy);

                    %%Step3b: Laplacian
                    [L_imxx, L_imxy] = gaussgradient(imx,obj.params.sigmagradient(1));
                    [L_imyx, L_imyy] = gaussgradient(imy,obj.params.sigmagradient(1));

                case 3
                    G = zeros(size(obj.img));
                    Diff_im = zeros(size(obj.img));
                    %Median filtering and diffusion in 2D. 
                    for i  = 1:size(obj.img,3)
                        G(:,:,i) = medfilt2(obj.img(:,:,i),obj.params.median_filter(1:2));
                        Diff_im(:,:,i) = diffusioncode(obj.img(:,:,i), obj.params.diffuse_iterations, 0.1429, obj.params.kappa1, obj.params.kappa2, obj.params.option);
                    end

                    %Rescale to [0,1]
                    Fim = mat2gray(Diff_im); clear Diff_im;
                    %Get gauss gradient. Using separable filters. sigmagradient now 3 component vector
                    [imx,imy]=gaussgradient2D_sep(Fim,obj.params.sigmagradient); clear Fim;
                    %Total magnitude. 
                    Mag_fim = hypot(imx,imy);

                    %Laplacian. 
                    [L_imxx, L_imxy] = gaussgradient2D_sep(imx,obj.params.sigmagradient);
                    [L_imyx, L_imyy] = gaussgradient2D_sep(imy,obj.params.sigmagradient);

            end

            Mag_Lim = L_imxx+L_imyy;  %Laplacian (trace of the 2D matrix)
            Mag_Lim=Mag_Lim.*(Mag_Lim>0);

            %%Step3c: Hessian Determniant
            Det_hessian=L_imxx.*L_imyy-L_imxy.*L_imyx;
            Det_hessian=Det_hessian.*(Det_hessian<0);
            %%Step4: Masking function using tanh on the Summed Derivatives from Step3
            X=obj.params.gamma-((obj.params.alpha*Mag_fim + obj.params.beta*Mag_Lim + obj.params.epsilon*abs(Det_hessian)) /obj.params.delta);
            Multi=0.5*(tanh(X)+1);% masking function
            obj.filtered = (double(G).*Multi); % masked image applied on smoothened image

        end
        
        
        function simple_thresholding(obj)
            % we refer to params.use_layer for thresholding. 
            obj.BW = obj.(obj.params.use_layer) >= obj.params.simple_threshold(obj.t);
            obj.bw2stats();

        end
        
        function hist_thresholding(obj)
            % we refer to params.use_layer for thresholding. 
            img = obj.(obj.params.use_layer);
            P = prctile(img(:),obj.params.percentile(obj.t));
            obj.BW = img >= P;
            obj.bw2stats();

        end
        
        function adaptive_thresholding(obj)
            % we refer to params.use_layer for thresholding. 
            img = obj.(obj.params.use_layer);

            %First find local maxima of filtered image. 
            BWmax = imextendedmax(img,obj.params.Hdepth);

            %Local max region props. 
            stats_max = regionprops(BWmax,img,'centroid','MeanIntensity','PixelIdxList');

            % Use the mean intensity of blobs to create an estimate of the best local
            % intensity threshold. 
            C = cat(1,stats_max.Centroid);
            ints = cat(1,stats_max.MeanIntensity);
            img_size = size(img);
            
            % 2d vs 3d. 
            switch length(img_size)
                
                case 2
                    reduction=8; % factor to reduce xy sampling. 
                    [xq,yq] = meshgrid(1:reduction:img_size(2),1:reduction:img_size(1));
                    F = scatteredInterpolant( C(:,1),C(:,2),ints, 'natural','nearest');
                    vq = F(xq,yq);
                    vq = imresize(vq,img_size);
                case 3
                    % This goes way faster by downsampling xy?
                    reduction=16; % factor to reduce xy sampling.                 
                    [xq,yq,zq] = meshgrid(1:reduction:img_size(2),1:reduction:img_size(1),1:img_size(3));
                    F = scatteredInterpolant( C(:,1),C(:,2),C(:,3),ints, 'natural','nearest');
                    % use parallel computing to run scattered interp. 
                    vq = zeros( size(xq) );

                    parfor i = 1:img_size(3)
                       vq(:,:,i) = F(xq(:,:,i),yq(:,:,i),zq(:,:,i)); 
                    end
                    %vq = F(xq,yq,zq);
                    % Now resize vq to actual image size. 
                    vq = imresize3(vq,img_size);
            end
            
            % Generate first estimate from vq and starting threshold. Uses
            % parameter: "local_threshold_factor" to obtain desired local
            % threshold. 
            obj.BW = img >= obj.params.local_threshold_factor.*vq;

            %obj.BW = obj.BW.*obj.mask;
            obj.bw2stats();
            
        end
        
        function close_filter(obj)
            
            %obj.stats = regionprops(obj.BW,'Centroid','Area','PixelList','PixelIdxList');

            %Now use imclose to remove gaps. Operate on individual cells so
            %that we minimize connecting cells that are close together. 
            se= strel('disk',obj.params.imclose_r,8);
            new_BW= zeros(size(obj.BW));

            dims=length(size(obj.BW));
            img_size=size(obj.BW);

            switch dims

                case 2

                    for c = 1:length(obj.stats)
                        p = index_2_vector( img_size, obj.stats(c).PixelIdxList );
                        X = p(:,1);
                        Y = p(:,2);

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
                        ind = sub2ind(size(obj.BW),Y,X);

                        %Fill in the large image. 
                        new_BW(ind) = 1;            

                    end
                    obj.BW = logical(new_BW);

                case 3
                    

                    for c = 1:length(obj.stats)
                        p = index_2_vector( img_size, obj.stats(c).PixelIdxList );

                        X = p(:,1);
                        Y = p(:,2);
                        Z = p(:,3);

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
                        ind = sub2ind(size(obj.BW),Y,X,Z);
                        %Fill in the large image. 
                        new_BW(ind) = 1;
                    end

                    obj.BW = logical(new_BW);


            end
            
            %Update stats. 
            obj.bw2stats();
        end
        
        
        function fill_holes(obj)
           
            if length(size(obj.BW))==2
                obj.BW = imfill(obj.BW,'holes');
                
            elseif length(size(obj.BW))==3

                pad_size = 60;
                Z = size(obj.BW,3);
                %Now fill in some holes. 
                for i = 1:Z
                    %Also fill holes laterally. Vertical holes can be
                    %problematic if we just to imfill with 3D stack

                    %first pad 
                    pad_plane = padarray(obj.BW(:,:,i),[pad_size,pad_size],'symmetric','both');
                    pad_plane_fill = imfill(pad_plane,'holes');

                    obj.BW(:,:,i) = pad_plane_fill(pad_size+1:end-pad_size, pad_size+1: end-pad_size);
                end              
            end            
            
            %Update stats. 
            obj.bw2stats();
            
        end
        
        function watershed(obj)
           
            %Connect regions, find objects that are are big enough to be doubles. 
            %Con=bwconncomp(logical(obj.BW), obj.params.nconnBW);
            %obj.stats = regionprops3(Con,'Centroid','Volume','SurfaceArea','VoxelIdxList');
            
            % we should already have 3D stats from any previous step!
            %obj.stats = obj.get_3D_stats(obj.BW);
            S = obj.sphericity_calculation';
            
            %Select the objects to watershed based off sphericity and
            %volume. 
            volumes = [obj.stats.Volume]';   
            liar_list = find( volumes > obj.params.WaterShedMaxVol & S <= obj.params.max_sphericity);
            img_dim= length(size(obj.BW));

            switch img_dim


                case 2

                    %Make bw of clusters only
                    bw_clust = zeros(size(obj.BW));
                    idx = cat(1,obj.stats(liar_list).PixelIdxList);
                    bw_clust(idx) = 1;

                    %distance transform
                    imgDist = -bwdist(~bw_clust);

                    %Smoothing
                    imgDist = medfilt2(imgDist,[3,3]);

                    %Get seeds 
                    mask = imextendedmin(imgDist, obj.params.h_min_depth, obj.params.h_min_conn); 

                    %Seems like smaller neighborhood works better?
                    imgDist = imimposemin(imgDist,mask);
                    imgDist(~bw_clust) = -Inf;

                    %Perform marked watershed
                    imgLabel = watershed(imgDist);

                    %Remove background again.
                    imgLabel(~bw_clust) = NaN;
                    stats_fused = obj.get_3D_stats(imgLabel);
                    %Skip first entry. 
                    stats_fused = stats_fused(2:end);
                    %ignore empty objects? 
                    sel = [stats_fused.Volume]==0;
                    stats_fused=stats_fused(~sel);

                    %Now get rid of original regions that were fused. 
                    keep = setdiff([1:length(obj.stats)],liar_list);
                    obj.stats = obj.stats(keep);

                    %Now append the fused nuclei obj.stats. 
                    obj.stats = [obj.stats;stats_fused];

                    %Now we need to adjust obj.BW... 
                    obj.BW = false(size(obj.BW));

                    %Add in new regions
                    new_idx = cat(1,obj.stats.PixelIdxList);
                    obj.BW(new_idx) = 1;

                case 3


                    %Distance transform scales. 
                    %scales = [obj.params.px_size(1),obj.params.px_size(2), obj.params.px_size(3)*obj.params.z_effect];
                    scales = [1,1,obj.params.z_effect];
                    %Counter of new objects found by watershedding
                    new_objects = 0;

                    %Empty fused object imglabel
                    fused_imgLabel = uint16(zeros(size(obj.BW))); 

                    %Loop over liars
                    for i = 1:length(liar_list)

                        %This liar idx
                        idx = liar_list(i);
                        display(['Watershedding object ',num2str(i), ' out of ', num2str(length(liar_list)) ])
                        %Sub_img containing this liar_volume.                         
                        Px = index_2_vector(size(obj.BW), obj.stats(idx).PixelIdxList); 
                        %Pixel ranges. 
                        min_range = min(Px,[],1);
                        max_range = max(Px,[],1);

                        x_range = min_range(1):max_range(1);
                        y_range = min_range(2):max_range(2);
                        z_range = min_range(3):max_range(3);

                        %Create a sub image of obj.BW. 
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
                        imgDist = medfilt3(imgDist,[5,5,1]);    

                        %Get seeds  %%ORIGINAL params were 0.7,6. With medfilt3 = 5,5,5. 
                        mask = imextendedmin(imgDist,obj.params.h_min_depth,obj.params.h_min_conn); %Seems like smaller neighborhood works better?
                        imgDist = imimposemin(imgDist, mask);
                        imgDist(~sub_bw) = -inf;

                        %Perform marked watershed
                        imgLabel = watershed(imgDist);
                        %Replace background again?
                        imgLabel(~sub_bw) = NaN;

                        %Putting objects into new compiled obj.BW. Figure out how many
                        %objects exist after watershedding. Infer from the values of
                        %imgLabel. imglabel=0 or 1 is background or envelope. 

                        %Number of objects
                        num_objects = double(max(imgLabel(:))-1 ); %Subtract 1 cause 0 and 1 are nothing. 

                        %Ignore non-object pixels. 
                        sel = imgLabel > 1;
                        imgLabel(~sel) = 0;

                        %Shift values by the on going counter. 
                        imgLabel(sel) = imgLabel(sel) + new_objects;

                        %Now place the imgLabel into the fused obj.BW img back where it
                        %came from. 
                        fused_imgLabel(y_range, x_range, z_range ) = fused_imgLabel(y_range, x_range, z_range ) + uint16(imgLabel);

                        %Adjust counter. 
                        new_objects = new_objects + num_objects;

                    end

                    %Keep non liar obj.stats
                    keep = setdiff([1:length(obj.stats)],liar_list);
                    obj.stats = obj.stats(keep); 

                    %Collect obj.stats on watershedded regions. 
                    stats_fused = obj.get_3D_stats(fused_imgLabel);
                    stats_fused = stats_fused(2:end); %Ignore first entry. 
                    %Ignore any empty areas. 
                    empty_list = [stats_fused.Volume]==0;
                    stats_fused = stats_fused(~empty_list);
                    %Append to obj.stats. 
                    obj.stats = [obj.stats; stats_fused];
                    %Now we need to adjust obj.BW... 
                    obj.BW = false(size(obj.BW));
                    %Add in new regions
                    new_idx = cat(1,obj.stats.PixelIdxList);
                    obj.BW(new_idx) = 1;

            end            
        end
        
        function merger(obj)
            
             %Centroids. 
            all_ctrs = cat(1,obj.stats.Centroid); 

            %Pairwise distance. 
            D = pdist2(all_ctrs,all_ctrs);

            %Select pair-wise distances less than threshold. 
            S = D < obj.params.merge_dist;

            %Get pairs. 
            [a,b] = find( tril( S ) );

            %Entries to get rid of. 
            bad_list = false(1,length(obj.stats));
            c=0;
            new_stats = struct('Centroid',[],'PixelList',[],'PixelIdxList',[],'Area',[]);
            
            %Merge algorithm. Loop over pairs, merge if result isn't too large.
            for m = 1:length(a)
                id_a = a(m);
                id_b = b(m);
                w_a = length(obj.stats(id_a).PixelIdxList);
                w_b = length(obj.stats(id_b).PixelIdxList);
                merge_vol = w_a + w_b;

                if merge_vol < obj.params.AbsMaxVol
                    bad_list(id_a) = 1;
                    bad_list(id_b) = 1;
                    new_list = [obj.stats(id_a).PixelIdxList; obj.stats(id_b).PixelIdxList];
                    new_centroid = w_a/merge_vol*obj.stats(id_a).Centroid + w_b/merge_vol*obj.stats(id_b).Centroid;
                    new_pixel_list = [obj.stats(id_a).PixelList; obj.stats(id_b).PixelList];
                    c=c+1;
                    new_stats(c).Centroid=new_centroid;
                    new_stats(c).PixelList = new_pixel_list;
                    new_stats(c).PixelIdxList = new_list;
                    new_stats(c).Area = length(new_pixel_list);
                end
            end

            %Now remove bad objects, and append new ones. 
            obj.stats = obj.stats(~bad_list);
            if isempty( new_stats(1).Centroid )
                disp('No merged objects')
            else
                obj.stats = [obj.stats; new_stats'];
            end
            obj.rebuild_BW();
            
        end
         
        
        function rg_filter(obj)
           
            
            rg = zeros(1,length(obj.stats));
            for i = 1:length(obj.stats)
               
                [y, x] = ind2sub(size(obj.BW),obj.stats(i).PixelIdxList);
                pos=[x,y];
                com=[mean(x),mean(y)];
                
                rg(i) = sqrt(1/length(y).*sum( sum( (pos-com).^2, 2)));
                
                
            end
            
            % Select Rg. 
            sel = rg <= obj.params.rg_threshold;
            
            obj.stats = obj.stats(sel);
            obj.rebuild_BW;
        end
        
        % Filtering cell/nuclei objects. 
        function filter_objects(obj)
            
            if isfield(obj.stats,'Area')
                volumes = [obj.stats.Area]';
            else 
                volumes = [obj.stats.Volume]';
            end

            sel = volumes >= obj.params.AbsMinVol & volumes <= obj.params.AbsMaxVol;
            obj.stats = obj.stats(sel);
            
            % Re-make BW. 
            obj.rebuild_BW;
        end
        
        
        % Apply the object mask to the filtered img. 
        function apply_mask(obj)
            obj.filtered = obj.filtered.*obj.mask;
            
        end
        
        % Local segmentation within masked objects. 
        function sub_region_segmentation(obj)
           
            obj.stats = struct('cell_id',[],'Area',[],'PixelIdxList',[], 'Centroid',[]);
            %counter. 
            c=0;
            
            % Erosion strel
            SE = strel('disk',10); %< why is this here???
            % Now perform a local segmentation per nuclei/cell. Define S/N threshold. 1.2? 
            n_cells = length(obj.msk_pxs);
            img_dims = length( size(obj.img) );
            for i = 1:n_cells
                
                % First make a mask, and erode as specified. 
                [sub_mask, ~, ranges] = px_list_2_mask( obj.msk_pxs{i}, size(obj.img));
                sub_mask = imerode(sub_mask, SE);
                
                switch img_dims
                    case 2                      
                        sub_img = obj.img(ranges.y,ranges.x);
                    case 3
                        % Sometimes segmented objects are just one z-plane. Let's
                        % skip those from processing.
                        if length(ranges.z) == 1
                            continue
                        end
                        sub_img = obj.img(ranges.y,ranges.x,ranges.z);
                end
                
                %the pixels. 
                vals= sub_img(sub_mask);

                % Clip prctile range
                int_range = prctile(vals,[obj.params.clip_prctile(1),obj.params.clip_prctile(2)]);
                vals_sel = vals( vals >= int_range(1) & vals <= int_range(2));

                % Calculate mean, std.     
                int_mean = mean(vals_sel);
                int_std  = std(vals_sel);

                % Skip if mean nucleus intensity is not above threshold: 
                if int_mean < obj.params.nucleus_int_threshold
                    continue
                end

                % Mask out non-cell pixels. 
                sub_img = sub_img.*sub_mask;

                % Threshold.
                thresh_val = int_std*obj.params.snr_ratio + int_mean;
                bw = logical(sub_img >= thresh_val);

                % Now generate regionprops for this nucleus. 
                sub_props = regionprops(bw, sub_img,'PixelIdxList','Area','Centroid');
                
                % Filter sizes. 
                areas = [sub_props.Area];
                sel_area = areas >= obj.params.foci_size_filter(1) & areas <= obj.params.foci_size_filter(2);
                sub_props = sub_props(sel_area);
                
                %loop over all sub-regions, add to blob list. Convert
                %pixelidxlist to large image coordinates. 
                for j = 1:length(sub_props)
                    c=c+1;

                    % cell id. 
                    sub_props(j).cell_id= i;
                    
                    % add background intensity. 
                    %sub_props(j).mean_bg = int_mean;
                    
                    % Calculate sum_intensity above background and SNR. 
                    %sum_int = sum( sub_img( sub_props(j).PixelIdxList ) - int_mean );
                    %snr = (mean( sub_img( sub_props(j).PixelIdxList ) ) - int_mean )/ int_std;
                    %sub_props(j).sum_int = sum_int;
                    %sub_props(j).snr     = snr;
                    
                    % Now convert index to large image. 
                    sub_props(j).PixelIdxList = sub_2_large_ind( sub_props(j).PixelIdxList, ranges, size(sub_img), size(obj.img));
                    
                    % convert centroid to large image units. 
                    ctr = sub_props(j).Centroid;
                    ctr(1) = ctr(1) + ranges.x(1) -1;
                    ctr(2) = ctr(2) + ranges.y(1) -1;
                    if img_dims == 3
                        ctr(3) = ctr(3) + ranges.z(1) -1;
                    end
                    sub_props(j).Centroid = ctr;

                    % now append to main structure. 
                    obj.stats(c)=sub_props(j);
                end
            end
            
            % rebuild BW from stats. 
            obj.rebuild_BW();
            
        end

        %% spot segment functions. 

        % pre-smooth the raw image. Not considered "filtered" yet. 
        function pre_smooth(obj)
            obj.img = imgaussfilt3(obj.img, obj.params.smooth_sigma);
        end

        % haven't tested...
        function subtract_local_background(obj)
        

            % Loop over cell nuclei in BW. 
            n = length(obj.msk_pxs);
            
            for i = 1:n
                
                bw = false(size(obj.BW));
                
                bw( obj.msk_pxs{i} ) = 1;
                
                mean_val = mean( obj.img(bw));
                %median_val = median( img(bw));
                
                obj.img(bw) = obj.img(bw)-mean_val;
                
            end
            
            % Set the background to zero?
            bg = ~obj.mask;
            obj.img(bg) = 0;
        end

        % LoG filtering. 
        function log_filter(obj)

            % Determine dimension. 
            switch length(size(obj.img))
                
                case 2
                    obj.filtered = -log_filter_2D(obj.img, obj.params.spot_sigma);

                case 3
                    obj.filtered = -log_filter_3D(obj.img, obj.params.spot_sigma);
            end

        end

        % thresholding spots. 
        function spot_simple_threshold(obj)
                    
            % Apply mask and threshold at once. 

            obj.BW = obj.filtered.*obj.mask >= obj.params.simple_threshold;

            %Stats on original image. Filter our single pixels. 
            obj.bw2stats();

            obj.assign_objects_2_cells()


        end
        
        %percentile-based local spot threshold. 
        function local_spot_threshold(obj)
                
            % Loop over all segmented mask objects. 
            N = length(obj.msk_pxs); % msks is a cell array,
            
            %Empty image to fill
            BW = zeros(size(obj.img));
            
            for i = 1:N
            
                %Get cell mask. Could be 2D or 3D. 
                switch length(size(obj.img))
                    
                    case 2
                        %2D index. 
                        I2d = obj.msk_pxs{i};
                        %Convert to 3D mask.
                        bw = false( size(obj.img,1), size(obj.img,2));
                        bw(I2d) = 1;
                        %Expand to proper z-dimensions. 
                        bw =repmat( bw, [1, 1, size(obj.img,3)] );
                        %Get index in 3D. 
                        this_cell = find(bw);
                        
                    case 3
                        
                        this_cell = obj.msk_pxs{i};
                        
                end
                
                %Get distribution of intensities
                int_dist_filter = obj.filtered( this_cell );
                int_dist_img    = obj.img( this_cell );
            
                %Threshold using filtered image
                thresh_val_filter = prctile(int_dist_filter,obj.params.local_thresh_percentile(2));
                
                %Threshold using original image
                thresh_val_og     = prctile(int_dist_img,obj.params.local_thresh_percentile(1));
                
                %Get pixels above percentile thresholds. 
                sel = int_dist_filter >= thresh_val_filter & int_dist_img >= thresh_val_og;
                
                
                %Add pixels to mask. 
                BW( this_cell(sel)) = 1;
            end
            
            %Collect Stats. 
            obj.BW=logical(BW);
            obj.bw2stats();

        end

        %Used for foci / spot detection and filtering.
        function S = sphericity_calculation(obj)
            
            
            %compactness https://mcib3d.frama.io/3d-suite-imagej/plugins/Analysis/Shape/3D-Shape/
            C = 36*pi*[obj.stats.Volume].^2 ./ [obj.stats.SurfaceArea].^3;
            %Sphericity
            S = C.^(1/3);

            vol=[obj.stats.Volume];
            sa  =[obj.stats.SurfaceArea];

            S = (pi^(1/3)).*(6.*vol).^(2/3) ./ sa;
        end

        %Used for foci / spot detection and filtering.
        function sphericity_filter(obj)
            %compactness https://mcib3d.frama.io/3d-suite-imagej/plugins/Analysis/Shape/3D-Shape/
            C = 36*pi*[obj.stats.Volume].^2 ./ [obj.stats.SurfaceArea].^3;
            %Sphericity
            S = C.^(1/3);
            
            sel = S >= obj.params.min_sphericity & S <= obj.params.max_sphericity;

            obj.stats=obj.stats(sel);            
            obj.rebuild_BW;
        end

        function add_mask(obj, msk_pxs)
            % If msk
            obj.msk_pxs = msk_pxs;
            % obj.msk_centroids = msk_centroids; don't need anymore

            % build a mask binary. Use size of BW because that can be
            % re-scaled. 
            obj.mask = false(size(obj.BW));
            obj.mask( cat(1, msk_pxs{:}) ) = 1;
            
        end

        % Resize image based on [y,x,z] dimensions
        function reverse_scaling(obj,img_scale)
            obj.BW = imbinarize(imresize3(double(obj.BW),img_scale));
            obj.bw2stats();
        end


        % Assign objects to cells using the mask data. Find the maximum
        % overlap between segmented object and cell mask. Sometimes with
        % image rescaling, the overlap isn't 100%. 
        function assign_objects_2_cells(obj)
           
            for i = 1:length(obj.stats)
                spot_pixels = obj.stats(i).PixelIdxList;
                U=zeros(1,length(obj.msk_pxs));
                for j = 1:length(obj.msk_pxs)
                    this_U = intersect(spot_pixels, obj.msk_pxs{j});
                    U(j) = length(this_U(:));
                end
                
                %through error if no overlap. 
                if sum(U)==0
                    warning('No overlap detected between segmented object and cell masks')
                     obj.stats(i).assignment=0;
                else

                    [~, obj.stats(i).assignment] = max(U);
                end
            end

        end

        % erode or dilate the masks. use negative values for erosion!
        function refine_mask(obj)
            
            % Erosion strel
            SE = strel('disk',abs(obj.params.mask_erosion_size));
            
            if obj.params.mask_erosion_size < 0
                obj.mask = imerode(obj.mask,SE);

            else
                obj.mask = imdilate(obj.mask,SE);
            end
            
            new_stats=regionprops(obj.mask,'Centroid','PixelIdxList');


            % Now compare intersection of PixelIdxlists. 
            og_msk = obj.msk_pxs;
            for i = 1:length(og_msk)
                U = zeros(1,length(new_stats));
                for j = 1:length(new_stats)
                    U(j) =length(intersect(og_msk{i},new_stats(j).PixelIdxList));
                end
                [~,assignment(i)] = max(U);  
            end
            
            stats = new_stats(assignment);
            obj.msk_pxs = list_2_cell_array( stats, 'PixelIdxList');
            
        end

        % spot fitting parent function. calls sub-functions to fit image
        % data. 
        % 
        % Params needed:
        %params.bg_mask
        %params.search_radius

        function spot_fitting(obj)
            

            switch obj.params.fit_type

                case 'arbitary'
                    % raw image, object stats, binary background mask, fit params. 
                    obj.spot_fits = fit_arbitrary(obj.img, obj.stats, obj.mask, obj.params );

                case 'gaussian'

            end

        end

        function filter_spot_snr(obj)
            sel = [obj.spot_fits.snr] >= obj.params.spot_snr_threshold;
            obj.spot_fits = obj.spot_fits(sel);

            % rebuild  BW. 
            obj.BW = false(size(obj.BW));
            obj.BW( cat(1,obj.spot_fits(:).PixelIdxList) ) = 1;
            obj.bw2stats()

        end

    end

            
    methods (Access=private)
    
        function rebuild_BW(obj)
           
            obj.BW = false(size(obj.img));
            if isfield(obj.stats,'PixelIdxList')
                idx = cat(1,obj.stats.PixelIdxList);
            elseif isfield(obj.stats,'VoxelIdxList')
                idx = cat(1,obj.stats.VoxelIdxList);
            end
            obj.BW(idx) = 1;
            
        end
        
        % need to update stats from bw operations. 
        function bw2stats(obj)
         
            % default to 3D stats for all seg_types. 
            obj.stats = obj.get_3D_stats(obj.BW);
 
        end
        
        % Useful function for standardizing the stats structure for 3D
        % regionprops. 
        function stats = get_3D_stats(obj,BW)
            stats = table2struct(regionprops3(BW, 'VoxelIdxList','Centroid','SurfaceArea','EigenVectors','EigenValues','Volume'));
            [stats.PixelIdxList] = stats.VoxelIdxList;
            stats = rmfield(stats,'VoxelIdxList');
        end



    end
    
end

        
        
        