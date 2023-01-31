% class for performing segmentation steps. 
classdef spot_segmenter < handle  
    
    properties (SetAccess = private)
        
        params
        
        step
        
        BW
        
        img
        
        filtered
        
        stats
        
        frame_obj
        
        mask
        
        t % frame.
    end
    
    
    
    methods (Access=public)
        
        % Constructer
        function obj = spot_segmenter( I, frame_obj, params) 
            
            obj.img = I;
            obj.mask= frame_obj.BW;
            obj.frame_obj = frame_obj;
            obj.params=params;
        end
        
        % pre smooth
        function pre_smooth(obj)
           
            obj.img = imgaussfilt3(obj.img, obj.params.smooth_sigma);
            
        end
        
        % Subtract local background. 
        function img = subtract_local_background(obj)
            
            % Get stats. 
            stats = regionprops(obj.mask,'PixelIdxList');

            % Loop over cell nuclei in BW. 
            n = length(stats);
            
            img=zeros(size(obj.mask));
            
            for i = 1:n

                bw = false(size(obj.mask));

                bw( stats(i).PixelIdxList ) = 1;

                mean_val = mean( obj.img(bw));
                %median_val = median( img(bw));

                img(bw) = img(bw)-mean_val;

            end

            % Set the background to zero?
            bg = ~obj.mask;
            obj.filtered(bg) = 0;
            
        end
            
        % Log filter with local background
        function log_filter_local_bg(obj)
           
            img = obj.subtract_local_background;
            
            switch length(size(img))
                
                case 2
                    
                    obj.filtered = -log_filter_2D(obj.img, obj.params.spot_sigma(1));
                    
                case 3
                    
                    obj.filtered = -log_filter_3D(obj.img, obj.params.spot_sigma);
                    
            end
            
            
        end
        
        % Log filter with local background
        function log_filter(obj)
                       
            switch length(size(obj.img))
                
                case 2
                    
                    obj.filtered = -log_filter_2D(obj.img, obj.params.spot_sigma(1));
                    
                case 3
                    
                    obj.filtered = -log_filter_3D(obj.img, obj.params.spot_sigma);
                    
            end
            
            
        end
        
        % Simple threshold of log image
        function simple_threshold(obj)
            
            obj.BW = obj.filtered >= obj.params.threshold;
            
            obj.stats = regionprops(logical(obj.BW), 'Centroid', 'PixelIdxList', 'Area');
            obj.assign_to_cells;
        end
        
        function mask_spots(obj)
            
            % Check if BW exists. 
            if isempty(obj.mask)
                error('No cell mask.')
            end

            % Mask out detected regions using cell mask. 
            if length(size(obj.mask))==2 & length(size(obj.filtered))==3
                mask = repmat(obj.mask,[1,1,size(obj.filtered,3)]);
            elseif length(size(obj.mask))==3 & length(size(obj.filtered))==3
                mask = obj.mask;
            else
                error('Incompatible img dims')
            end
            
            % Mask out regions. 
            obj.BW = logical(obj.BW .* mask);
            obj.stats = regionprops(obj.BW,'Centroid','PixelIdxList','Area');  
            obj.assign_to_cells;
        end
        
        function filter_objects(obj)
            
            volumes = [obj.stats.Area]';
            
            sel = volumes >= obj.params.AbsMinVol & volumes <= obj.params.AbsMaxVol;
            obj.stats = obj.stats(sel);
            
            % Re-make BW. 
            obj.rebuild_BW;
        end
        
        function filter_sizes(obj)
            
            areas = [obj.stats.Area];
            
            sel = areas >= obj.params.sizes(1) & areas <= obj.params.sizes(2);
            obj.stats = obj.stats(sel);
            obj.rebuild_BW;
        end
    end
    
    methods (Access=private)
    
        function rebuild_BW(obj)
           
            obj.BW = false(size(obj.BW));
            idx = cat(1,obj.stats.PixelIdxList);
            obj.BW(idx) = 1;
            
        end
        
        function assign_to_cells(obj)
           
            % assign to cells. 
            spot_centroids = cat(1, obj.stats.Centroid);
            spot_seg_dim = length( obj.stats(1).Centroid);
            
            cell_centroids = cat(1,obj.frame_obj.centroids{:});
            cell_seg_dim = size(cell_centroids,2);
            
            if spot_seg_dim == cell_seg_dim
                D = pdist2(spot_centroids, cell_centroids);
            elseif spot_seg_dim ==2 && cell_seg_dim==3
                cell_centroids = [cell_centroids, h/2*ones(size(cell_centroids,1),1)];
                D = pdist2(spot_centroids, cell_centroids);
            elseif spot_seg_dim ==3 && cell_seg_dim==2
                D = pdist2(spot_centroids(:,1:2), cell_centroids);
            end

            [~,assignment] = min(D,[],2);
            A = num2cell(assignment);
            [obj.stats.assignment] = A{:};
            
        end
        
    end
    
end

        
        
        