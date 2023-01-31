% Calculate mean yap signal (could be arbitrary signal based on image).
% Uses the roi2poly tool. Maybe faster?
%Now computes cyto to nuclear ratio by taking a mean of a thin ring outside
%of the contour

%***Updated to use structures for calculations. 

%***2-27-18 Now we will search through the stack to find the best plane(s) for
%calculating N/C signals. We will find which plane gives the maximum
%nuclear intensity for the pre-calculated mask. Then take that plane +/- a
%few other planes to do MIP and nuc/cyto calculations

function obj = nuc_cyto_2D(obj, params, step, force_frames)


step   = default_step( step );

debug = 0;
if(debug)
    figure(7)
    clf(7)
    %set(gcf,'color','w')    
    color_vec  = hsv(length(obj.tracks));
    color_vec  = color_vec(randperm(size(color_vec,1)),:);
end


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
    %Create empty image to fill up
    msk_img = zeros(size_y,size_x,size_z);
    sig_img = zeros(size_y,size_x,size_z);
    
    %Load the Frame_obj file associated with this image
    if ~exist( seg_files{t})
        continue
    end
    
    load(seg_files{t});
    
    %Get the bio-formats image index corresponding to this z-stack:
    msk_img = get_stack(reader,t,params.seg_channel);
    sig_img = get_stack(reader,t,params.sig_channel);
   
    
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

        % If mask is 3D, then project into 2D. Get center from Z-location
        % of 3D segmented cell.
        mask_dims = length(size(frame_obj.(channel_str).BW));
        if mask_dims == 3
            BW = zeros(size_y,size_x,size_z);
            px = frame_obj.(channel_str).PixelIdxList{j};
            BW(px)=1;
            % Max project.
            BW = max(BW,[],3);
            %Get pixel coordinates. 
            [px_Y, px_X, px_Z] = ind2sub([size_y,size_x,size_z],px);

            %Create mask. 
            roi_bw_nuc = logical(BW);
            
            idx = round(frame_obj.(channel_str).centroids{j}(3));
            
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

            % Figure out optimal z-section. 
            % Get sub_img of nucleus. 
            msk_range_y = [min(px_Y):max(px_Y)];
            msk_range_x = [min(px_X):max(px_X)];
            sub_nuc = msk_img(msk_range_y,msk_range_x,:);
            sub_mask = roi_bw_nuc(msk_range_y,msk_range_x,:);
            sub_mask = repmat( sub_mask,[1,1,Z]);
            sub_nuc(~sub_mask) = NaN;

            %Find mean signal in each plane. 
            plane_vals = reshape(  mean( mean( sub_nuc,1,'omitnan'),2,'omitnan'),[1,Z]);

            %Find index of max value <-i.e. the 'best' plane
            [~,idx] = max(plane_vals);

        end
        

        %% Make a range of planes to use based on params.z_range
        z_planes = idx-params.z_range: idx+params.z_range;
        %Need to adjust z-planes if they are outside of img
        if(z_planes(end) > Z)
            sel = z_planes <= Z;
            z_planes = z_planes(sel);
        end
        if( z_planes(1) < 1 )
            sel = z_planes >= 1;
            z_planes = z_planes(sel);
        end
        
        
        %% Analyze the signal channel. 
        %Look at sub-region of image containing this cell. Need to exapnd
        %to include outer boundary. 
        x_min = max(1,min(px_X)-params.out_boundary);
        x_max = min(size_x,max(px_X)+params.out_boundary);
        y_min = max(1,min(px_Y)-params.out_boundary);
        y_max = min(size_y,max(px_Y)+params.out_boundary);
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
        
        %Use intensity of nuclear region to minimize nucleolar regions.
        if step.nuc_thresh

            %Which img do we use? 
            if strcmp(params.nuc_thresh_channel,'self')
                nuc_plane=sig_img(:,:,idx);
            elseif stcmp(params.nuc_thresh_channel,'seg_channel')
                nuc_plane = msk_img(:,:,idx);
            end
            
            %Set non-nuclear region to NaN. 
            nuc_plane(~roi_bw_nuc) = NaN;
            T = prctile(nuc_plane(:),params.nuc_prctile);
            nuc_mask = nuc_plane >= T;
            int_mask = nuc_mask(y_range,x_range);
        else
            int_mask = roi_bw_nuc(y_range,x_range);
        end
        
        %The clean mask - represents boundaries of nucleus (filled). 
        clean_nuc_mask = roi_bw_nuc(y_range,x_range);                
        
        %Collect data. 
        data(j) = analyze_roi(sub_img, clean_nuc_mask, int_mask, j, z_planes, params );
        
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

%% ROI analysis function.

function data = analyze_roi(img, clean_mask, int_mask, ID, z_planes , params )

%% Nucleus. 

%First check if clean_mask and int_mask are equal.
if isequal(clean_mask,int_mask)
    
    %Erode mask a by in-boundary pixels.First pad with zeros. << forgot why
    %padding?
    mask_pd = padarray(clean_mask,[1,1],0);
    
    %If we need buffer away from segmentation edge. 
    if( params.in_buffer>0)
        se = strel('disk',params.in_buffer,8);
        mask_outter = imerode(mask_pd,se);
    else
        mask_outter = mask_pd;
    end
    
    se = strel('disk',params.in_boundary,8);
    mask_inner = imerode(mask_outter,se);
    %unpad.
    mask_inner = mask_inner(2:end-1,2:end-1);
    mask_outter = mask_outter(2:end-1,2:end-1);
    %Use some mask differences to query the ROIs. 
    roi_inner = logical( mask_outter - mask_inner);
else
    
    roi_inner = int_mask;
    
    
end

%% Cytoplasm. 
    %Now get cyto mask.
    
    %Buffer as needed
    if(params.out_buffer > 0)
        se = strel('disk',params.out_buffer,8);
        mask_inner = imdilate(clean_mask,se);
    else
        mask_inner = clean_mask;
    end
    se = strel('disk',params.out_boundary,8);
    mask_outer = imdilate(mask_inner,se);

    roi_outer = logical( mask_outer - mask_inner);
    

    %% Calculate the mean vals of nuclear and cytoplasmic regions

    %Make all masks the correct depth in 3D.
    roi_inner = repmat(roi_inner,[1,1,length(z_planes)]);
    roi_outer   = repmat(roi_outer,[1,1,length(z_planes)]);

    %Get mean/median nuclear and cytoplasm fluorescence 
    data.nuc_mean   = mean( img( roi_inner ),'omitnan' );
    data.nuc_med    = median( img( roi_inner ),'omitnan' );
    data.cyto_mean  = mean( img( roi_outer   ),'omitnan' );   
    data.cyto_med   = median( img(roi_outer),'omitnan');
    data.cell_id    = ID;
    data.local_rho  = NaN;
    data.nuc_area   = sum( clean_mask(:) );
    if(isnan(data.nuc_mean) || isnan(data.cyto_mean) || isnan(data.cyto_med) || isnan(data.nuc_med))
        disp('warning: nan found')
        disp(['cell: ',num2str(ID)]);
    end


%% DEBUGGING code for plotting contours. 
debug=0;
if debug 
    c = contourc(double(max(roi_inner,[],3)),[0.5,0.5]);
    C = C2xyz(c);
    newFigure(18)
    imshow3D(img);    
        hold on

    for i = 1:length(C)
        plot(C{i}(:,1),C{i}(:,2),'linewidth',2)
    end
end



end

%% Generate an empty structure for keeping track of nuclear and cytoplasmic signals
function data = gen_data_struct( N )

%Fields
data(N) = struct('nuc_mean',[],'cyto_mean',[],'local_rho',[],'cell_id',[],'nuc_med',[],'cyto_med',[],'nuc_area',[]);

end



%% defaults. 
function step = default_step( step )

%List of all default parameters. 
dstep.nuc_thresh=0;

S  = fieldnames( dstep );

for i = 1:length(S)
    
    %Check if this field exists. 
    if ~isfield(step,S{i})
        step.(S{i}) = dstep.(S{i});
        %Output this default was used. 
        disp(['Using default ',S{i},' with value: ',num2str(dstep.(S{i}))]);
    end
end

end

