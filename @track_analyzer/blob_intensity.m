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
    % seg_type: 'foci', 'cells', etc. 
    % sig_channel:
    % dilation_buffer:
    % params.background_distance:
    % params.background_thickness:

% Optional
    % max_p = True -> use max-projection image. must already be made. 

if ~isfield(params,'max_p')
    params.max_p = false;
end

% Channel str.
channel_str = ['seg_channel_',pad(num2str(params.seg_channel),2,'left','0')];

T = obj.exp_info.t_frames;

if nargin < 4
    frames = 1:T;
else
    frames = force_frames;
end

%check if we use max projection. 
if params.max_p
    %Generate reader. FOR NOW, assume we are looking in series 1. 
    reader = bfGetReader(obj.exp_info.max_p_img);
    series = 1;
else
    %Generate reader. FOR NOW, assume we are looking in series 1. 
    reader = bfGetReader(obj.exp_info.img_file);
    series = 1;
end

%Get the image size of this series. 
size_x = reader.getSizeX;
size_y = reader.getSizeY;
size_z = reader.getSizeZ;

% make a big BW mask. 
if size_z > 1
    BW = false(size_y,size_x,size_z);
else
    BW = false(size_y,size_x);
end

%Get segmentation file names. 
seg_files = obj.get_frame_files;

disp(['Intensity analysis using channel ',channel_str])

for t = frames
    display(['Analyzing frame: ',num2str(t)])
    tic

    %Load the Frame_obj file associated with this image
    if ~exist( seg_files{t})
        continue
    end
    
    load(seg_files{t});
    
    %Get the bio-formats image index corresponding to this z-stack:
    msk_img = get_stack(reader,t,params.seg_channel);
    sig_img = get_stack(reader,t,params.sig_channel);
   
    % Get data to perform calculations. 

    seg_data = frame_obj.(channel_str).(params.seg_type);
    % Get segmentation stats. rebuild a bw. 
    stats = seg_data.stats;
    bw = obj.rebuild_BW(cat(1,stats.PixelIdxList));
    
    % dilate all masks by 'dilation_buffer' pixels. Use this to
    % prevent overlap of local background estimation between distinct objects. 
    se = strel('disk',params.dilation_buffer,8);
    bw_all_masks_dilated = imdilate(bw,se);


    %Clear out bad values (stitching can create 0 valued pixels at edge). 
    zSEL = sig_img == 0;
    sig_img(zSEL) = NaN;
    zSEL = msk_img == 0;
    msk_img(zSEL) = NaN;
   
    %Loop over detected objects
    if isfield(seg_data.stats,'PixelIdxList')
        n_objs = length(seg_data.stats);
    else
        disp(['No objects found in frame: ', num2str(t)])
        continue
    end
            
    %Create empty structure for nuclear / cytoplasm calculations
    data = gen_data_struct( n_objs );

    if usejava('desktop')
        f = waitbar(0, 'Starting');
    end


    %Loop over segmented objects. 
    for j = 1:n_objs
        
        % wait bar. 
        if usejava('desktop')
            waitbar(j/n_objs, f, sprintf('Progress: %d / %d', j, n_objs));
            pause(0.001);        
        end

        % 3D vs 2D 
        mask_dims = length(size(bw));
        if mask_dims == 3
            
            %Force BW to be zero again (saves time this way)
            roi_bw_nuc = false(size(BW));
            px = seg_data.stats(j).PixelIdxList;
            roi_bw_nuc(px)=true;
            %Get pixel coordinates. 
            [px_Y, px_X, px_Z] = ind2sub([size_y,size_x,size_z],px);
                        
        % Estimate z-plane from max intensity.
        elseif mask_dims == 2

            %Get a 2D mask of this cell....
            roi_bw_nuc = false(size(BW));
            px = seg_data.stats(j).PixelIdxList;
            roi_bw_nuc(px) = true;
            
            %Get pixel coordinates. 
            [px_Y,px_X] = ind2sub([size_y,size_x],px);

        end
        

        %% Make a range of planes to use based on params.z_range
        if size_z == 1
            z_planes = 1;
        else
            z_planes = 1:size_z;
        end
        
        
        %% Analyze the signal channel. 
        %Look at sub-region of image containing this object. Need to exapnd
        %to include outer boundary.  
        out_boundary = params.background_distance + params.background_thickness;
        x_min = max(1,min(px_X)-out_boundary);
        x_max = min(size_x,max(px_X)+out_boundary);
        y_min = max(1,min(px_Y)-out_boundary);
        y_max = min(size_y,max(px_Y)+out_boundary);
                   
        %New [x,y] range. 
        x_range = [x_min:x_max];
        y_range = [y_min:y_max];

        if mask_dims == 3
            z_min = max(1,min(px_Z)-out_boundary);
            z_max = min(size_z,max(px_Z)+out_boundary);
            z_range = [z_min:z_max];
        end

        sub_img = sig_img(y_range,x_range,z_range);
        
        %Tricky part: we need to ignore pixels belonging to other nuclei! 
            %All PixelIdx NOT belonging to cell j. 
            sel = [1:n_objs] ~= j;
            % 3D case.
            if size_z > 1
                ind = cat(1,seg_data.stats(sel).PixelIdxList);
                not_this_cell = false(size_y,size_x,size_z);
                not_this_cell(ind) = true;
                nan_mask = not_this_cell(y_range,x_range,z_range);
            else
                ind = cat(1,seg_data.stats(sel).PixelIdxList);
                not_this_cell = false(size_y,size_x);
                not_this_cell(ind) = 1;
                %3D mask.
                not_this_cell = repmat(not_this_cell,[1,1,length(z_planes)]);
                nan_mask = not_this_cell(y_range,x_range,:);
            end
            
            %multiply sub_stack by the mask 
            sub_img(nan_mask)= NaN;
            
        %Now get sub_img using sub_stack. Still a img volume though...
        %sub_img = double(sub_stack(y_range,x_range,:));
        sub_mask = roi_bw_nuc(y_range,x_range,z_range);
        sub_bw_all_masks_dilated = bw_all_masks_dilated(y_range,x_range,z_range);
        
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
        data(j).local_background = mean(bg_vals);
%         
%         figure(1); clf
%         imshow3D(sub_background_mask)
%         figure(2); clf
%         imshow3D(sub_bw_all_masks_dilated)
% 
%         pause

       % figure(1); imshow3D(sub_img)
       % figure(2); imshow3D(not_this_cell(y_range,x_range,:))
    end

    %Now save the frame_obj. Append new variables to existing structure.
    fnames=fieldnames(data);

    % check if the segmented channel has data saved already
    seg_channel = ['seg_channel_',pad(num2str(params.sig_channel),2,'left','0')];
    if isfield(frame_obj, seg_channel)
        og_data = frame_obj.(seg_channel).(params.seg_type);
        
        for i = 1:length(fnames)
            [og_data.stats.(fnames{i})] = data.(fnames{i});
        end
        frame_obj.(['seg_channel_',pad(num2str(params.sig_channel),2,'left','0')]).(params.seg_type) = og_data;
    else
        % Make a new field in frame_obj for this channel. 
        frame_obj.(seg_channel).(params.seg_type).stats = data;
    end
    save(seg_files{t},'frame_obj','-append')
    
    
    if usejava('desktop')    
        delete(f);
    end

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
data(N) = struct('sum_int',[],'mean_int',[],'cell_id',[],'local_background',[]);

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
function debugge()
c = contourc(double(max(sub_background_mask,[],3)),[0.5,0.5])
C=C2xyz(c)
imshow3D(img)
hold on
for i = 1:length(C)
    plot(C{i}(:,1),C{i}(:,2),'linewidth',2)
end

end