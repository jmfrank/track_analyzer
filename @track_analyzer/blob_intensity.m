
% Calculate mean signal inside a segmented obj. 

% Use dilation to calculate local background values. 

function obj = blob_intensity(obj, params, force_frames)


debug = 0;
if(debug)
    figure(7)
    clf(7)
    %set(gcf,'color','w')    
    color_vec  = hsv(length(obj.tracks));
    color_vec  = color_vec(randperm(size(color_vec,1)),:);
end

% parameters needed: 
    % seg_channel:
    % sig_channel:
    % dilation_buffer:
    % params.background_distance:
    % params.background_thickness:

% Channel str.
channel_str = ['seg_channel_',pad(num2str(params.seg_channel),2,'left','0')];

Z = obj.exp_info.z_planes;
T = obj.exp_info.t_frames;

if nargin < 4
    frames = 1:T;
else
    frames = force_frames;
end

%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = bfGetReader(obj.exp_info.img_file);
series = 1;

%Get the image size of this series. 
size_x = reader.getSizeX;
size_y = reader.getSizeY;
size_z = reader.getSizeZ;

%Get segmentation file names. 
seg_files = obj.get_frame_files;

for t = frames
    display(['Analyzing mean signal of frame: ',num2str(t)])
    tic

    %Load the Frame_obj file associated with this image
    if ~exist( seg_files{t})
        continue
    end
    
    load(seg_files{t});
    
    %Get the bio-formats image index corresponding to this z-stack:
    msk_img = get_stack(reader,t,params.seg_channel);
    sig_img = get_stack(reader,t,params.sig_channel);
   
    % get bw masks of whole image. 
    bw = frame_obj.(channel_str).BW;
    
    % dilate all masks by 'dilation_buffer' pixels. Use this to
    % prevent overlap of local background estimation between distinct objects. 
    se = strel('disk',params.dilation_buffer,8);
    bw_all_masks_dilated = imdilate(bw,se);


    %Clear out bad values (stitching can create 0 valued pixels at edge). 
    zSEL = sig_img == 0;
    sig_img(zSEL) = NaN;
    zSEL = msk_img == 0;
    msk_img(zSEL) = NaN;
   
    %Loop over detected cells
    if isfield(frame_obj.(channel_str),'PixelIdxList')
        n_cells = length(frame_obj.(channel_str).PixelIdxList);
    else
        disp(['No cells in frame: ', num2str(t)])
        continue
    end
            
    %Create empty structure for nuclear / cytoplasm calculations
    data = gen_data_struct( n_cells );

    %Loop over cells. 
    f = waitbar(0, 'Starting');

    %Loop over cells. 
    for j = 1:n_cells
        
        % wait bar. 
        waitbar(j/n_cells, f, sprintf('Progress: %d %%', floor(j/n_cells*100)));
        pause(0.001);        

        % this is old code---need to re-write. 
        mask_dims = length(size(frame_obj.(channel_str).BW));
        if mask_dims == 3
            
            BW = zeros(size_y,size_x,size_z);
            px = frame_obj.(channel_str).PixelIdxList{j};
            BW(px)=1;

            %Get pixel coordinates. 
            [px_Y, px_X, px_Z] = ind2sub([size_y,size_x,size_z],px);

            %Create mask. 
            roi_bw_nuc = logical(BW);
                        
        % Estimate z-plane from max intensity.
        elseif mask_dims == 2

            %Get a 2D mask of this cell....
            BW = zeros(size_y,size_x);
            px = frame_obj.(channel_str).PixelIdxList{j};
            BW(px) = 1;
            
            %Get pixel coordinates. 
            [px_Y,px_X] = ind2sub([size_y,size_x],px);

            %Create mask. 
            roi_bw_nuc = logical(BW);

        end
        

        %% Make a range of planes to use based on params.z_range
        if size_z == 1
            z_planes = 1;
        else
            z_planes = 1:size_z;
        end
        
        
        %% Analyze the signal channel. 
        %Look at sub-region of image containing this cell. Need to exapnd
        %to include outer boundary. We don't want to look at pixels above
        %and below blobs for background. 
        out_boundary = params.background_distance + params.background_thickness;
        x_min = max(1,min(px_X)-out_boundary);
        x_max = min(size_x,max(px_X)+out_boundary);
        y_min = max(1,min(px_Y)-out_boundary);
        y_max = min(size_y,max(px_Y)+out_boundary);
        
        %New [x,y] range. 
        x_range = [x_min:x_max];
        y_range = [y_min:y_max];
               
        sub_stack = sig_img(:,:,z_planes);
        
        %Tricky part: we need to ignore pixels belonging to other nuclei! 
            %All PixelIdx NOT belonging to cell j. 
            sel = [1:n_cells] ~= j;
            % 3D case.
            if size_z > 1
                ind = cat(1,frame_obj.(channel_str).PixelIdxList{sel});
                not_this_cell = false(size_y,size_x,size_z);
                not_this_cell(ind) = true;
                not_this_cell = not_this_cell(:,:,z_planes);
            else
                ind = cat(1,frame_obj.(channel_str).PixelIdxList{sel});
                not_this_cell = false(size_y,size_x);
                not_this_cell(ind) = 1;
                %3D mask.
                not_this_cell = repmat(not_this_cell,[1,1,length(z_planes)]);
            end
            
            %Using not_this_cell as mask, replace all pixels with NaN. 
            sub_stack( not_this_cell ) = NaN;
            
        %Now get sub_img using sub_stack. Still a img volume though...
        sub_img = double(sub_stack(y_range,x_range,:));
        sub_mask = roi_bw_nuc(y_range,x_range,:);
        sub_bw_all_masks_dilated = bw_all_masks_dilated(y_range,x_range,:);
        
        % first dilate to background distance. 
        se_inner = strel('disk',params.background_distance,8);
        sub_mask_background_inner = imdilate(sub_mask,se_inner);
        se_outer = strel('disk',params.background_distance + params.background_thickness,8);
        sub_mask_background_outer = imdilate(sub_mask,se_outer);
        % Difference in masks. 
        sub_background_mask = logical(sub_mask_background_outer - sub_mask_background_inner);
        % Need to multiply by inverse of sub_bw_masks_dilated. 
        sub_background_mask = logical(sub_background_mask .* (~sub_bw_all_masks_dilated));
        % Get sum, mean intensity of masked region. 
        vals = sub_img(sub_mask);
        bg_vals = sub_img(sub_background_mask);
        data(j).sum_int = sum(vals);
        data(j).mean_int = mean(vals);
        data(j).cell_id = j;
        data(j).nuc_area = sum(sub_mask(:));
        data(j).local_background = mean(bg_vals);
        
%         figure(1); clf
%         imshow3D(sub_background_mask)
%         figure(2); clf
%         imshow3D(sub_bw_all_masks_dilated)
% 
%         pause
    end
    close(f);

    %Now save the frame_obj. Saving to a channel specific field.
    frame_obj.(['channel_',pad(num2str(params.sig_channel),2,'left','0')]) = data;
    save(seg_files{t},'frame_obj','-append')

    if(debug)
        pause
        if(i < length(obj.img_files));clf(7);
        end
    end
    toc
end

if(debug);hold off;end

end

%% Generate an empty structure for keeping track of nuclear and cytoplasmic signals
function data = gen_data_struct( N )

%Fields
data(N) = struct('sum_int',[],'mean_int',[],'cell_id',[],'nuc_area',[],'local_background',[]);

end
%% DEBUGGING code for plotting contours. 
% debug=0;
% if debug 
%     c = contourc(double(max(roi_inner,[],3)),[0.5,0.5]);
%     C = C2xyz(c);
%     newFigure(18)
%     imshow3D(img);    
%         hold on
% 
%     for i = 1:length(C)
%         plot(C{i}(:,1),C{i}(:,2),'linewidth',2)
%     end
% end
% 
% 
% 
% end


