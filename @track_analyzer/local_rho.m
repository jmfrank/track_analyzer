%Calculate local cell density. performed using binary operations based on
%cell centroid pixels and circular masks centered at each cell. 

function [obj] = local_rho(obj,params)
    
%Generate empty mask the same size as image. 
BW = false(obj.exp_info.img_size);

%Create a structuring element for the search area. 
SE = strel('disk',params.search_radius,8);

%Now loop over time, load frame data, assess density, save. 
F = obj.get_frame_files;
for i = 1:length(F)
    
    disp_str = [];
    %display(['Analyzing frame #',num2str(i)])   

    load( F{i} );
    
    %Figure out data type. Append if existing structure. 
    if isfield(frame_obj,'channel_01')
        data = frame_obj.channel_01;
        out_field = 'channel_01';
    else
        data = struct('local_rho',[]);
        out_field = 'data';
    end
    
    %Centroids. 
    centroids = round( cat(1,frame_obj.centroids{:}));
    
    %Loop over cells. 
    n_cells = length(frame_obj.centroids);
    for c = 1:n_cells
        
        %Start with mask for cell of interest. 
        mask = BW;
        ctr = centroids(c,:);
        mask(ctr(2),ctr(1)) = 1;
        
        %Dilate. 
        mask = imdilate(mask,SE);
        
        %Now see how many centroids fit in dilated mask. 
        sel = [1:n_cells] ~= c;
        px = sub2ind(size(BW),centroids(sel,2),centroids(sel,1));
        K = BW;
        K(px) = 1;
        
        %Density. Cells per square-pixels. 
        data(c).local_rho = sum( K(mask) ) ./ sum( mask );
    end
    
    %Append data. 
    frame_obj.(out_field) = data;
    
    
end


end



