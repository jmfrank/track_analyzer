% Generate a mask from pixel list. 

function [mask, large_mask, ranges] = px_list_2_mask( p , img_size )

% dimensionality. 
switch length(img_size)
    
    case 2
        
        % Check what format p is in. If it's just indices or coordinates. 
        if isvector(p)

            [pY,pX] = ind2sub(img_size, p);
            p = [pX,pY];

        end

        %Pixel ranges. 
        min_range = min(p,[],1);
        max_range = max(p,[],1);

        x_range = min_range(1):max_range(1);
        y_range = min_range(2):max_range(2);

        %Create a sub image of BW. 
        mask = false(length(y_range),length(x_range));

        %Add in the ones from this blob. 
        X = p(:,1) - x_range(1) + 1;
        Y = p(:,2) - y_range(1) + 1;

        ind = sub2ind(size(mask),Y,X);
        mask(ind) = 1;

        % Large mask of original image. 
        large_mask = false(img_size);
        idx = sub2ind(img_size, p(:,2),p(:,1));
        large_mask(idx) = 1;

        % ranges. 
        ranges.x=x_range;
        ranges.y=y_range;
    
    case 3

        % Check what format p is in. If it's just indices or coordinates. 
        if isvector(p)

            [pY,pX,pZ] = ind2sub(img_size, p);
            p = [X,Y,Z];

        end


        %Pixel ranges. 
        min_range = min(p,[],1);
        max_range = max(p,[],1);

        x_range = min_range(1):max_range(1);
        y_range = min_range(2):max_range(2);
        z_range = min_range(3):max_range(3);

        %Create a sub image of BW. 
        mask = false(length(y_range),length(x_range),length(z_range));

        %Add in the ones from this blob. 
        X = p(:,1) - x_range(1) + 1;
        Y = p(:,2) - y_range(1) + 1;
        Z = p(:,3) - z_range(1) + 1;

        ind = sub2ind(size(mask),Y,X,Z);
        mask(ind) = 1;

        % Large mask of original image. 
        large_mask = false(img_size);
        idx = sub2ind(img_size, p(:,2),p(:,1),p(:,3));
        large_mask(idx) = 1;

        % ranges. 
        ranges.x=x_range;
        ranges.y=y_range;
        ranges.z=z_range;
    
end