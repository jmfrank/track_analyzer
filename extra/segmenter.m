% class for performing segmentation steps. 
classdef segmenter < handle  
    
    properties (SetAccess = private)
        
        params
        
        step
        
        BW
        
        img
        
        filtered
        
        stats
        
        image_bits
        
        t % frame.
    end
    
    
    
    methods (Access=public)
        
        % Constructer
        function obj = segmenter( I, image_bits, params, t ) 
            
            obj.img = I;
            obj.image_bits= image_bits;
            obj.params=params;
            obj.t=t;
        end
        
        
        % Image processing functions. 
        function mean_filter(obj)
            
            % Check bit-depth. 
            max_bits = max(obj.img(:));
            if max_bits <= 256
                image_bits=8;
            elseif max_bits <= 4096
                image_bits=12;
            elseif max_bits <= 65536
                image_bits=16;
            end
            
            if image_bits==8
                if ~strcmp(class(obj.img),'uint8')
                    obj.img = uint8(obj.img);
                end
                tmp = adaptthresh(obj.img,obj.params.MeanFilterSensitivity,'NeighborhoodSize',obj.params.MeanFilterNeighborhood);
            elseif image_bits==12

                %Rescale to 16 bit? 
                obj.img = uint16(obj.img* (2^16-1)/(2^12-1));
                tmp = adaptthresh(obj.img,obj.params.MeanFilterSensitivity,'NeighborhoodSize',obj.params.MeanFilterNeighborhood);

            elseif image_bits==16
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
           
            obj.BW = obj.filtered >= obj.params.simple_threshold;
            [obj.BW, obj.stats] = filter_objects(obj.BW,obj.params.AbsMinVol);

        end
        
        function hist_thresholding(obj)
           
            P = prctile(obj.filtered(:),obj.params.percentile(obj.t));
            
            obj.BW = obj.filtered >= P;
            [obj.BW, obj.stats] = filter_objects(obj.BW,obj.params.AbsMinVol);

        end
        
        function adaptive_thresholding(obj)
           
            p_val = obj.params.percentile(obj.t);

            %Threshold starting point. 
            obj.params.thresh_start = prctile(obj.filtered(:),p_val);
                if obj.params.thresh_start == 0
                    warning('Threshold too low');
                    %Calculate the first non-zero percentile and set as
                    %starting percentile. 
                    ps = prctile(obj.filtered(:),[p_val:0.2:100]);
                    nonzero = ps > 0;
                    obj.params.thresh_start = nonzero(1);

                end

            %Incremental change in threshold. 
            %params.increment = params.thresh_start*0.01;
            %First smooth image. 
            I_sm = imgaussfilt(obj.img,obj.params.I_sm_sigma);
            %Perform iterative thresholding. 
            obj.BW = iterative_thresholding(I_sm, obj.filtered, obj.params );
            [obj.BW, obj.stats] = filter_objects(obj.BW,obj.params.AbsMinVol);

        end
        
        function close_filter(obj)
            
            obj.stats = regionprops(obj.BW,'Centroid','Area','PixelList','PixelIdxList');

            %Now use imclose to remove gaps. Operate on individual cells so
            %that we minimize connecting cells that are close together. 
            se= strel('disk',obj.params.imclose_r,8);
            new_BW= zeros(size(obj.BW));

            dims=length(size(obj.BW));

            switch dims

                case 2

                    for c = 1:length(obj.stats)
                        X = obj.stats(c).PixelList(:,1);
                        Y = obj.stats(c).PixelList(:,2);

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
                        X = obj.stats(c).PixelList(:,1);
                        Y = obj.stats(c).PixelList(:,2);
                        Z = obj.stats(c).PixelList(:,3);

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
        end
        
        function watershed(obj)
           
            %Connect regions, find objects that are are big enough to be doubles. 
            Con=bwconncomp(logical(obj.BW), obj.params.nconnBW);
            obj.stats = regionprops(Con,'Centroid','Area','PixelList','PixelIdxList');

            %New volume list
            volumes = [obj.stats.Area]';   
            liar_list = find( volumes > obj.params.WaterShedMaxVol);
            disp('Starting watershed')


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
                    stats_fused = regionprops(imgLabel,'Centroid','Area','PixelList','PixelIdxList');

                    %Skip first entry. 
                    stats_fused = stats_fused(2:end);

                    %Now get rid of original regions that were fused. 
                    keep = setdiff([1:length(obj.stats)],liar_list);
                    obj.stats = obj.stats(keep);

                    %Now append the fused nuclei obj.stats. 
                    obj.stats = [obj.stats;stats_fused];

                    %Now we need to adjust obj.BW... 
                    obj.BW = zeros(size(obj.BW));

                    %Add in new regions
                    new_idx = cat(1,obj.stats.PixelIdxList);
                    obj.BW(new_idx) = 1;

                    disp('Finished watershed')

                case 3


                    %Distance transform scales. 
                    scales = [obj.params.px_size(1),obj.params.px_size(2), obj.params.px_size(3)*obj.params.z_effect];
                    scales = [1,1,obj.params.z_effect];
                    %Counter of new objects found by watershedding
                    new_objects = 0;

                    %Empty fused object imglabel
                    fused_imgLabel = uint16(zeros(size(obj.BW))); 

                    %Loop over liars
                    for i = 1:length(liar_list)

                        %This liar idx
                        idx = liar_list(i);

                        %Sub_img containing this liar_volume. 
                        Px = obj.stats( idx ).PixelList;

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
                        imgDist = medfilt3(imgDist,[3,3,3]);    

                        %Get seeds  %%ORIGINAL params were 0.7,6. With medfilt3 = 5,5,5. 
                        mask = imextendedmin(imgDist,obj.params.h_min_depth,obj.params.h_min_conn); %Seems like smaller neighborhood works better?
                        imgDist = imimposemin(imgDist,mask);
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
                    stats_fused = regionprops(fused_imgLabel,'Centroid','Area','PixelList','PixelIdxList');
                    stats_fused = stats_fused(2:end); %Ignore first entry. 
                    %Ignore any empty areas. 
                    empty_list = [stats_fused.Area]==0;
                    stats_fused = stats_fused(~empty_list);
                    %Append to obj.stats. 
                    obj.stats = [obj.stats; stats_fused];
                    %Now we need to adjust obj.BW... 
                    obj.BW = zeros(size(obj.BW));
                    %Add in new regions
                    new_idx = cat(1,obj.stats.PixelIdxList);
                    obj.BW(new_idx) = 1;

                    disp('Finished watershed')

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

            %Merge algorithm. Loop over pairs, merge if result isn't too large.
            for m = 1:length(a)
                id_a = a(m);
                id_b = b(m);
                w_a = length(obj.stats(id_a).PixelIdxList);
                w_b = length(obj.stats(id_b).PixelIdxList);
                merge_vol = w_a + w_b;

                if merg_vol < AbsMaxVol
                    bad_list(id_a) = 1;
                    bad_list(id_b) = 1;
                    new_list = [obj.stats(id_a).PixelIdxList; obj.stats(id_b).PixelIdxList];
                    new_centroid = w_a/merge_vol*obj.stats(id_a).Centroid + w_b/merge_vol*obj.stats(id_b).Centroid;
                    new_pixel_list = [obj.stats(id_a).PixelList; obj.stats(id_b).PixelList];
                    c=c+1;
                    new_stats(c).Centroid=new_centroid;
                    new_stats(c).PixelList = new_pixel_list;
                    new_stats(c).PixelIdxList = new_list;
                end
            end

            %Now remove bad objects, and append new ones. 
            obj.stats = obj.stats(~bad_list);
            obj.stats = [obj.stats; new_stats];
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
        
        function filter_objects(obj)
            
            volumes = [obj.stats.Area]';
            
            sel = volumes >= obj.params.AbsMinVol & volumes <= obj.params.AbsMaxVol;
            obj.stats = obj.stats(sel);
            
            % Re-make BW. 
            obj.rebuild_BW;
        end
    end
    
    methods (Access=private)
    
        function rebuild_BW(obj)
           
            obj.BW = false(size(obj.BW));
            idx = cat(1,obj.stats.PixelIdxList);
            obj.BW(idx) = 1;
            
        end
        
    end
    
end

        
        
        