% Calculate mean yap signal (could be arbitrary signal based on image).
% Uses the roi2poly tool. Maybe faster?
%Now computes cyto to nuclear ratio by taking a mean of a thin ring outside
%of the contour

%***Updated to use structures for calculations. 

%***2-27-18 Now we will search through the stack to find the best plane(s) for
%calculating N/C signals. We will find which plane gives the maximum
%nuclear intensity for the pre-calculated mask. Then take that plane +/- a
%few other planes to do MIP and nuc/cyto calculations

function obj = nuc_cyto_signal_analysis(obj, params, step, force_frame)


debug = 0;
if(debug)
    figure(7)
    clf(7)
    %set(gcf,'color','w')    
    color_vec  = hsv(length(obj.tracks));
    color_vec  = color_vec(randperm(size(color_vec,1)),:);
end



Z = obj.exp_info.z_planes;
T = obj.exp_info.t_frames;

if nargin < 4
    frames = 1:T;
else
    frames = force_frame;
end

%UPDATED: USING BIOFORMATS READER
%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = bfGetReader(obj.exp_info.img_file);
series = 1;

%Get the image size of this series. 
size_x = reader.getSizeX;
size_y = reader.getSizeY;

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
    for i = 1:Z
        this_plane = reader.getIndex(i-1,params.seg_channel-1,t-1)+1;
        msk_img(:,:,i) = bfGetPlane(reader,this_plane);
        
        this_plane = reader.getIndex(i-1,params.signal_channel-1,t-1)+1;
        sig_img(:,:,i) = bfGetPlane(reader,this_plane);
        
    end
    
    %Loop over detected cells
    try
        n_cells = length(frame_obj.refined_cells);
        seg_type = 'refined';
    catch
        if(~isfield(frame_obj,'PixelIdxList'))
            continue
        end
        n_cells = length(frame_obj.PixelIdxList);
        seg_type = 'Px';
    end
    
    %Create empty structure for nuclear / cytoplasm calculations
    data = gen_data_struct( n_cells );
        
    %Loop over cells. 
    for j = 1:n_cells
        
        %Get a 2D mask of this cell, depending on segmentation format....
        switch seg_type
            
            %DEPRECATED.....
            case 'refined'
                %Contour
                x     = frame_obj.refined_cells{j}(:,1);
                y     = frame_obj.refined_cells{j}(:,2);
                %Get mask of this cell
                roi_bw_nuc = poly2mask(x,y,size_y,size_x);  
            %SHOULD ALWAYS use Px. 
            case 'Px' 
                BW = zeros(size_y,size_x);
                px = frame_obj.PixelIdxList{j};
                BW(px) = 1;
                %Get pixel coordinates. 
                [px_Y,px_X] = ind2sub([size_y,size_x],px);
                
                %Create mask. 
                roi_bw_nuc = logical(BW);
        end
        
        %% Figure out optimal z-section. 
        %Get sub_img of nucleus. 
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
        %Make a range of planes to use based on params.z_range
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
        
        %Sub-sample image to get proper z-range. 
        sub_stack = sig_img(:,:,z_planes);
        
        %Tricky part: we need to ignore pixels belonging to other nuclei! 
            %All PixelIdx NOT belonging to cell j. 
            sel = [1:n_cells] ~= j;
            ind = cat(1,frame_obj.PixelIdxList{sel});
            not_this_cell = false(size_y,size_x);
            not_this_cell(ind) = 1;
            %3D mask.
            not_this_cell = repmat(not_this_cell,[1,1,length(z_planes)]);       
            %Using not_this_cell as mask, replace all pixels with NaN. 
            sub_stack( not_this_cell ) = NaN;
            
        %Now get sub_img using sub_stack. Still a img volume though...
        sub_img = double(sub_stack(y_range,x_range,:));
        
        %Use intensity of DNA to minimize nucleolar regions.
        if step.nuc_thresh
            nuc_plane = msk_img(:,:,idx);
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
    
    %Now save the frame_obj
    frame_obj.data = data;
    save(seg_files{t},'frame_obj','-append')

    if(debug);pause;if(i < length(obj.img_files));clf(7);end;end
    toc
end

       
if(debug);hold off;end;

obj.nuc_cyto_calc = true;

%Now 
end


%% Generate an empty structure for keeping track of nuclear and cytoplasmic signals
function data = gen_data_struct( N )

%Fields
data(N) = struct('nuc_mean',[],'cyto_mean',[],'local_rho',[],'cell_id',[],'nuc_med',[],'cyto_med',[],'nuc_area',[]);

end




