%Calculate local cell density. performed using binary operations based on
%cell centroid pixels and circular masks centered at each cell. 

function [obj] = local_rho(obj,params)
    
%Generate empty mask the same size as image. 
BW = false(obj.exp_info.img_size);

%Create a structuring element for the search area. 
%SE = strel('disk',params.search_radius,8);

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
    %All centroids in px coordinates.
    px = sub2ind(size(BW),centroids(:,2),centroids(:,1));
    K = BW;
    K(px) = 1;
    
    %Loop over cells. 
    n_cells = length(frame_obj.centroids);
    for c = 1:n_cells
        
        %Start with mask for cell of interest. 
        mask = BW;
        ctr = centroids(c,:);
                
        %Range. 
        x_min= max(1,ctr(1)-params.search_radius);
        y_min = max(1,ctr(2)-params.search_radius);
        x_max = min(obj.exp_info.img_size(2),ctr(1)+params.search_radius);
        y_max = min(obj.exp_info.img_size(1),ctr(2)+params.search_radius);
        
        sub_bw = BW(y_min:y_max,x_min:x_max);
        
        %Location of all pixels. 
        [Y,X]=find( ~sub_bw );
        
        %dist. 
        new_ctr = ctr - [x_min,y_min];
        D = ((Y-new_ctr(2)).^2 + (X-new_ctr(1)).^2).^0.5;
        
        sel = D <= params.search_radius;
        
        sub_bw(sel) = 1;
        
        %mask(y_min:y_max,x_min:x_max)=sub_bw;
        
        sub_K = K(y_min:y_max,x_min:x_max);
        
        
        %Dilate. 
        %mask = imdilate(mask,SE);
        
        %Density. Cells per square-pixels. 
        n_cells = sum( sub_K(sub_bw) ) - 1; %Subtract self. 
        area_px = sum( sub_bw(:) );
        data(c).local_rho =  n_cells / area_px ;
    end
    
    %Append data. 
    frame_obj.(out_field) = data;
    
    %Save. 
    save(F{i},'frame_obj');
    
    if exist('disp_str','var'); clearString(disp_str);end
    disp_str=['Frame ',num2str(i)];
    disp(disp_str)  
end

%Add local rho value. 
obj.exp_info.search_radius=params.search_radius;
obj.save;

end



