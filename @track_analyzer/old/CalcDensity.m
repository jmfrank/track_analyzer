%Updating this to use the updated structures in obj.

%For each time point, estimate the local density of cells for each contour.
%We need to remove cells that are near the boundaries of the image.
%Specifically, params.density_dist is the search distance, so cells should
%be at least that far away from any edge. Need to open a sample image and
%get image size to make proper calculation. 

function [obj, total_density] = CalcDensity(obj,params)
 


%Open sample image
IMG = tiffread(obj.exp_info.img_file,1);
height = IMG.height;
width = IMG.width;


%Looping over time points ( length of obj.cells )
skip = 0; %Counter for skipping cells
for i = 1:length(obj.cells)
    %display(['Analyzing frame #',num2str(i)])   

    %get x,y position of each cell
    dat = cell2mat(cellfun(@(x) mean(x), obj.cells{i},'UniformOutput',false)');
    D = pdist2(dat,dat);

    %Now looping over 1st dimension of D (corresponds to the cells found at
    %frame 'i'
    for j = 1:size(D,1)
        
        %This cell's center position
        this_cell = dat(j,:);
        %Check if it's sufficiently far from the edge. Pretty sure first
        %column is 'x' and second column is 'y'. CHECK THIS
        if( this_cell(1) < params.density_dist || this_cell(1) > width - params.density_dist ...
                || this_cell(2) < params.density_dist || this_cell(2) > height - params.density_dist )
            
            %Cell is too close to boundary
            skip  = skip + 1;
            neighbors(i,j) = NaN; %Don't count these cells
            
        else
            %Cell is sufficiently far from image boundary

            %Find the number of neighbors within the distance metric, subtract one for
            %counting self distance
            neighbors(i,j) = sum( D(j,:) <= params.density_dist) - 1;
            %Now find the nearest neighbor
            sel = D(j,:) > 0;
            nn(i,j)        = min( D(j,sel));
        end 
    end
    
    %Calculate the total density of each frame. Output is in cells / um^2? 
    total_density(i) = size(dat,1) ./ (obj.exp_info.pixel_size^2*width*height);
end

%Now go through the tracking data and add in the neighbor metric to the
%nuc_cyto_data structure
for i = 1:length(obj.tracks)
    for j = 1:size(obj.tracks{i},1)
        %Frame
        t = obj.tracks{i}(j,1);
        %Cell ID
        c = obj.tracks{i}(j,2);
        %Find the ID in the nuc_cyto_data structure
        data_id = obj.tracks{i}(j,5);
        %Add density to the local_rho field as cells/um^2
        obj.nuc_cyto_data(data_id).local_rho = neighbors(t,c) / pi / (obj.exp_info.pixel_size*params.density_dist)^2;
    end
end


end




