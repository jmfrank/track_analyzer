%Local thresholding on each cell. Using a histogram binning method to find
%pixels that are above a certain percentile of the intensity distribution
%found in each cell. This helps account for cell-2-cell expression
%variation. 

%% 1-29-23 update. Need to send frame_obj info specifically from the nuclei-segmentation channel. 
% frame_obj has two fields, PixelidxList (cell array, as from output of bw
% improp calculations), and centroids. 
function [stats] = local_threshold_cells( img, img_filter, frame_obj, indices, params)

%Loop over cells. 
n_cells = length(indices);

%Dimensions
[w, l, h] = size(img);

%Empty image to fill
BW = zeros(size(img));

for i = 1:n_cells

    %Get cell mask. Could be 2D or 3D. 
    switch params.seg_dims
        
        case 2
            %2D index. 
            I2d = frame_obj.PixelIdxList{indices(i)};
            %Convert to 3D mask.
            bw = false( size(img,1), size(img,2));
            bw(I2d) = 1;
            %Expand to proper z-dimensions. 
            bw =repmat( bw, [1, 1, size(img,3)] );
            %Get index in 3D. 
            this_cell = find(bw);
            
        case 3
            
            this_cell = frame_obj.PixelIdxList{indices(i)};
            
            
    end
    
    %Get distribution of intensities
    int_dist_filter = img_filter( this_cell );
    int_dist_img    = img( this_cell );

    %Threshold using filtered image
    thresh_val_filter = prctile(int_dist_filter,params.local_thresh_percentile(2));
    
    %Threshold using original image
    thresh_val_og     = prctile(int_dist_img,params.local_thresh_percentile(1));
    
    %Get pixels above percentile thresholds. 
    sel = int_dist_filter >= thresh_val_filter & int_dist_img >= thresh_val_og;
    
    %If we have pixels, continue. 
    if(sum(sel)==0)
    %    continue
    end
    
    %Add pixels to mask. 
    BW( this_cell(sel)) = 1;
end

%Collect Stats. 
stats = table2struct(regionprops3(logical(BW), 'VoxelIdxList','Centroid','SurfaceArea','EigenVectors','EigenValues','Volume'));
[stats.PixelIdxList] = stats.VoxelIdxList;
stats = rmfield(stats,'VoxelIdxList');

%Assignment. Spot must reside within cell. but which one? Take a single
%pixel and find intersection to cell pixels? Or make a bwlabel?
for i = 1:length(stats)

    this_spot_pixel = stats(i).PixelIdxList(1);
    for j = 1:length(frame_obj.PixelIdxList)
        U = intersect(this_spot_pixel, frame_obj.PixelIdxList{j});
        if ~isempty(U)
            stats(i).assignment = j;
        end
    end
end


end


%%


% for c = 1:n_cells
%     
%     ctr = frame_obj.(channel_str).centroids{c};
%     
%     text(ctr(1),ctr(2),num2str(c),'color','w')
% end


%     %% Plotting and other functionalities to ignore for now. 
%     function plot_stuff
%     
%     %Build images
%     x_range = min(Xn)-1:max(Xn)+1;
%     y_range = min(Yn)-1:max(Yn)+1;
%     z_range = min(Zn)-1:max(Zn)+1;    
%     x_range = x_range >=1 & x_range <= l;
%     y_range = y_range >=1 & y_range <= w;
%     z_range = z_range >=1 & z_range <= h;
%     sub_img = img(y_range,x_range,z_range);
%     sub_img_filt = img_filter(y_range,x_range,z_range);
%     
%     %Binarize
%     sub_img_mask = sub_img >= thresh_val_og;
%     sub_img_mask_filt = sub_img_filt >= thresh_val_filter;
%     
%     %Loop over each bright pixel and evaluate separately. 
%     for p = 1:length(Yn)
%         Y=Yn(p);
%         X=Xn(p);
%         Z=Zn(p);
%         %Neighboring indices
%         neighbors = [Y-1,X,Z;
%                      Y+1,X,Z;
%                      Y,X-1,Z;
%                      Y,X+1,Z;
%                      Y,X,Z-1;
%                      Y,X,Z+1];
%         %Make sure they didn't go out of bounds
%         sel_n = neighbors(:,1) >= 1 & neighbors(:,1) <= w & neighbors(:,2) >=1 & neighbors(:,2)<=l & neighbors(:,3)>=1 & neighbors(:,3)<=h;
%         neighbors = neighbors(sel_n,:);
%         
%         %Get intensities
%         ind_neighbors = sub2ind(size(img),neighbors(:,1),neighbors(:,2),neighbors(:,3));
%         int_neighbors = img(ind_neighbors);
% 
%         %Is at least one neighbor 1std above mean val? 
%     end
%                      
%     %Plotting
%     [Yn,Xn,Zn] = ind2sub(size(img),this_cell);
%     x_range = min(Xn):max(Xn);
%     y_range = min(Yn):max(Yn);
%     z_range = min(Zn):max(Zn);
%     
%     sub_img_raw = img(y_range,x_range,z_range);   
%     sub_img_filt= img_filter(y_range,x_range,z_range);
%     
%     %Bright pixels
%     bright_pixels = [Yn(sel),Xn(sel),Zn(sel)];
%     %Shift to sub_img
%     bright_pixels = bright_pixels - [min(Yn),min(Xn),min(Zn)]+1;
%     %Put as ind
%     IND = sub2ind(size(sub_img_raw),bright_pixels(:,1),bright_pixels(:,2),bright_pixels(:,3));
%     sub_img_mask = zeros(size(sub_img_raw));
%     sub_img_mask(IND) = 1;
%     
%     global_max_val = max(sub_img_raw(:));
%     bg_filter = -log_filter_3D(sub_img_raw,[15,15,10]);
%     rescale_filter = sub_img_filt - bg_filter;
%     max_val = max(rescale_filter(:));
%     min_val = min(rescale_filter(:));
%     
%     rescale_filter = (rescale_filter - min_val) ./ (max_val - min_val)*global_max_val;
%     sub_img_mask = global_max_val*sub_img_mask;
%     C = horzcat(sub_img_raw,rescale_filter,sub_img_mask);
%     
%     figure(2); clf;
%     imshow3D(C);
%     pause
%     
%     end
    




%