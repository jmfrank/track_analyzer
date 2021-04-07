% Calculate mean yap signal (could be arbitrary signal based on image).
% Uses the roi2poly tool. Maybe faster?
%Now computes cyto to nuclear ratio by taking a mean of a thin ring outside
%of the contour

%***Updated to use structures for calculations. 

%***2-27-18 Now we will search through the stack to find the best plane(s) for
%calculating N/C signals. We will find which plane gives the maximum
%nuclear intensity for the pre-calculated mask. Then take that plane +/- a
%few other planes to do MIP and nuc/cyto calculations

function obj = nuc_cyto_3D(obj, params, step, force_frames)

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
    msk_img = zeros(size_y,size_x,Z);
    sig_img = zeros(size_y,size_x,Z);
    
    %Load the Frame_obj file associated with this image
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
    for j = 1:n_cells
        % wwait bar. 
        waitbar(j/n_cells, f, sprintf('Progress: %d %%', floor(j/n_cells*100)));
        pause(0.1);
        
        %Get a 3D mask of this cell....
        BW = zeros(size_y,size_x,size_z);
        px = frame_obj.(channel_str).PixelIdxList{j};
        BW(px) = 1;
        
        %Get pixel coordinates. 
        [px_Y, px_X, px_Z] = ind2sub([size_y,size_x, size_z],px);

        %% Analyze the signal channel. 
        %Look at sub-region of image containing this cell. Need to exapnd
        %to include outer boundary. 
        x_min = max(1,min(px_X)-params.out_boundary(1));
        x_max = min(size_x,max(px_X)+params.out_boundary(1));
        y_min = max(1,min(px_Y)-params.out_boundary(2));
        y_max = min(size_y,max(px_Y)+params.out_boundary(2));
        z_min = max(1,min(px_Z)-params.out_boundary(3));
        z_max = min(size_z,max(px_Z)+params.out_boundary(3));
        
        %New [x,y] range. 
        x_range = [x_min:x_max];
        y_range = [y_min:y_max];
        z_range = [z_min:z_max];
        
        %Create mask. 
        roi_bw_nuc = logical(BW);
        sub_mask = roi_bw_nuc(y_range,x_range, z_range);
        
        % Temporary store a copy of the signal image for masking. 
        this_img = sig_img;
        
        %Tricky part: we need to ignore pixels belonging to other nuclei! 
            %All PixelIdx NOT belonging to cell j. 
            sel = [1:n_cells] ~= j;
            ind = cat(1,frame_obj.(channel_str).PixelIdxList{sel});
            not_this_cell = false(size_y,size_x,size_z);
            not_this_cell(ind) = 1;
            %Using not_this_cell as mask, replace all pixels with NaN. 
            this_img( not_this_cell ) = NaN;
            
        %Now get sub_img using sub_stack. Still a img volume though...
        sub_img = double(this_img(y_range, x_range, z_range));
        
        %Use intensity of nuclear region to minimize nucleolar regions.
        if step.nuc_thresh

            %Which img do we use? 
            if strcmp(params.nuc_thresh_channel,'self')
                nuc_plane=sub_img;
            elseif stcmp(params.nuc_thresh_channel,'seg_channel')
                nuc_plane = msk_img(y_range, x_range, z_range);
            else 
                error('No nuclear threshold channel defined');
            end
            
            %Set non-nuclear region to NaN. 
            nuc_plane(~sub_mask) = NaN;
            T = prctile(nuc_plane(:), params.nuc_prctile);
            int_mask = nuc_plane >= T;
        else
            int_mask = sub_mask;
        end
                
        %Collect data. 
        data(j) = analyze_roi(sub_img, sub_mask, int_mask, j, params );
        %kkk=1
    end
    close(f);
    
    %Now save the frame_obj. Saving to a channel specific field.
    frame_obj.(['channel_',pad(num2str(params.sig_channel),2,'left','0')]) = data;
    save(seg_files{t},'frame_obj','-append')

    if(debug);pause;if(i < length(obj.img_files));clf(7);end;end
    toc
end

if(debug);hold off;end

%obj.nuc_cyto_calc = true;

%Now 
end

%% Analyze NC. 

function data = analyze_roi(img, clean_mask, int_mask, ID, params )

%% Nucleus. 

%First check if clean_mask and int_mask are equal.
if isequal(clean_mask,int_mask)
    
    %Erode mask a by in-boundary pixels.First pad with zeros. << forgot why
    %padding?
    mask_pd = padarray(clean_mask,[1,1,1],0);
    
    %If we need buffer away from segmentation edge. 
    if( params.in_buffer>0)
        se = strel('sphere',params.in_buffer);
        mask_outter = imerode(mask_pd,se);
    else
        mask_outter = mask_pd;
    end
    
    se = strel('sphere',params.in_boundary);
    mask_inner = imerode(mask_outter,se);
    %unpad.
    mask_inner = mask_inner(2:end-1,2:end-1,2:end-1);
    mask_outter = mask_outter(2:end-1,2:end-1,2:end-1);
    %Use some mask differences to query the ROIs. 
    roi_inner = logical( mask_outter - mask_inner);

else
    
    roi_inner = int_mask;
    
    
end

%% Cytoplasm. 
    %Now get cyto mask.
    
    %Buffer as needed
    if(params.out_buffer > 0)
        se = strel('sphere',params.out_buffer);
        mask_inner = imdilate(clean_mask,se);
    else
        mask_inner = clean_mask;
    end
    se = strel('sphere',params.out_boundary(1));
    mask_outer = imdilate(mask_inner,se);

    roi_outer = logical( mask_outer - mask_inner);
    

    %% Calculate the mean vals of nuclear and cytoplasmic regions

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




