%This is a function for creating a frap trace. Specific to nascent
%transcription frap experiments. 

function obj = frap_spot_tracking_v2(obj, params)

debug = 0;

%Initiate track. 
frap_track = [];

%% Check if there's a pre-bleach image. 
if( isfield( obj.exp_info,'pre_bleach') )

    %We know the first frame must contain a nascent transcript. Load first
    %frame file. 
    reader = bfGetReader(obj.exp_info.pre_bleach);
    pre_bleach_T = params.pre_bleach_T;
    X = reader.getSizeX;
    Y = reader.getSizeY;
    Z = reader.getSizeZ;
    T_pre = reader.getSizeT;
    IMG_pre = zeros(Y,X,Z,T);
    for t = 1:T_pre
        IMG_pre(:,:,:,t) = get_stack(reader,t,1);
    end
    reader.close;

else
    T_pre = 0;
end

%% Now load img file. 
reader = bfGetReader(obj.exp_info.img_file);
X = reader.getSizeX;
Y = reader.getSizeY;
Z = reader.getSizeZ;
T = reader.getSizeT;
IMG_all = zeros(Y,X,Z,T+T_pre);
if T_pre > 0
    IMG_all(:,:,:,1:T_pre) = IMG_pre;
end

for t = T_pre+1:T
    IMG_all(:,:,:,t) = get_stack(reader,t,1);
end
reader.close;

%Total number of frames. 
n_frames = size(IMG_all,4);
%Number of pre-bleach frames. 
pre_bleach_T = params.pre_bleach_T;

%% Deal with some ROI things...get roi. Force to use the ellipse roi. 
roi_types = cellstr(strvcat(obj.exp_info.ROI.type));
roi = obj.exp_info.ROI( cellfun(@(x) strcmp(x,'Ellipse'),roi_types));

%if there's a rectangle ROI, mask out the image outside of rectangle. 
[yes,rec_idx] = ismember('Rectangle',roi_types);
if(yes)
    plane = logical(zeros(Y,X));
    plane(obj.exp_info.ROI(rec_idx).PixelIdx) = 1;
    cols = find(any(plane,1));
    rows = find(any(plane,2));
else
    rows = 1:Y;
    cols = 1:X;
end

%Shift roi to proper position. 
roi.x_pos = roi.x_pos - cols(1);
roi.y_pos = roi.y_pos - rows(1);

%% Loop over time for pre-bleach image. Collect info. 
fits_all = gen_fit_struct(pre_bleach_T);
IMG = zeros(length(rows),length(cols),n_frames);
for t = 1:pre_bleach_T
    
    %This stack. 
    stack = IMG_all(:,:,:,t);
    
    %Fit frame. 
    fit = get_fits(stack(rows,cols,:), params);      
    pos = cat(1,fit.pos);
    pos = pos(:,[2,1]);
    
    %Dist to roi center. 
    D = sum( (pos - [roi.x_pos,roi.y_pos]).^2, 2).^0.5;
    %Find fit_id that's closest. 
    [~, min_idx] = min(D);
    %Create fits_all object with closest spot as first entry. 
    fits_all(t) = fit(min_idx);
    %add first frame to track (t=1, fit id, position)
    frap_track(t,:) = [t,t,pos(min_idx,:)];
    
    %Build up max-p image. 
    IMG(:,:,t) = max(stack(rows,cols,:),[],3);  

end

%% We know the spot is totally gone in first post-bleach frame, so measure the intensity where the spot was. 
stack = IMG_all(:,:,:,pre_bleach_T+1);
%Project. 
mip = max(stack(rows,cols,:),[],3);
%Append to image. 
IMG(:,:,end+1) = mip;

%Query region where spot was in pre-bleach image. Make a 'fake' stats structure. 
stats.Centroid = [fits_all(end).pos(2),fits_all(end).pos(1)];
stats.assignment = 1;
switch params.fit_type
    case '2D'
        
        fit = fit_this_frame_centroid_2D(mip,stats,params);
    case '3D'
        fit = fit_this_frame_centroid(stack,stats,params);
end

%Add this fit to track. 
fits_all(end+1) = fit;
frap_track(end+1,:) = [pre_bleach_T+1,length(fits_all),fit.pos(2),fit.pos(1)];

%Now run through the rest of the time points and track. 

count = length(fits_all); %Count the total number of fits to this point. 
xyzt = []; %Build up for tracking. 
for t = pre_bleach_T+2:n_frames
    
    stack = IMG_all(:,:,:,t);
    %Append to IMG. 
    IMG(:,:,end+1)=max(stack(rows,cols,:),[],3);
    
    %Get fits
    fit = get_fits(stack(rows,cols,:),params);
    
    if(isempty(fit))
        continue
    end
    
    %Filter out fit that's not bright.
    sel = cat(1,fit.sum_int) >= params.sum_int_thresh;
    fit = fit(sel);
    
    if(isempty(fit))
        continue
    else
              
       %Position vector. comes as y,x, (z)
        pos = cat(1,fit.pos);
        %Switch to x,y. Ignoring z. 
        pos = pos(:,[2,1]);

        %ID's for results?
        ids = [1:size(pos,1)] + count;
        %Append fits to xyzt. [position, ids, time (frames)]
        xyzt = [xyzt; pos, ids', t.*ones(size(pos,1),1)];
        %Append fits_all structure
        fits_all = [fits_all,fit];
        %Update counter. 
        count = length(fits_all);
    end
end

%Now try tracking. 
tracks = track(xyzt,params.max_disp,params);

%% Reorganize data. Tracks are: pos, ids, time. We want spot_tracks as: time, ids, pos. 
    ids= unique(tracks(:,end));
    spot_tracks = cell( length(ids),1);
    %New columns
    d=length(tracks(1,:));
    new_columns = [d-1,d-2,1:size(pos,2)];
    for i = 1:length(ids)
        sel = tracks(:,end)==ids(i);
        spot_tracks{i} = tracks(sel,new_columns);
    end

    % Now try to link one of the spot_tracks to the data. Use minimum distance
    %to first point in track.
    %first_pos = cell2mat(cellfun(@(x) x(1,3:4),spot_tracks,'uniformoutput',0));
    %D = sum((fits_all(1).pos([2,1]) - first_pos).^2,2).^0.5;
    %[~,min_idx] =min(D);
    %check which time point this track starts. 
    %t_start = spot_tracks{min_idx}(1,1);

%Pre-pend spot-tracks with the bleach track. 
spot_tracks = [frap_track;spot_tracks];

%Open viewer for connecting tracks manually. Receive track ids, and fits. 
track_ids = imshow2D_time_lapse(IMG,spot_tracks);

%Now find when frames are missing and interpolate image. 
cat_track = cell2mat( spot_tracks( track_ids ) );

%Now we have to look if there's any overlaps in times. 
all_t = unique( cat_track(:,1));
corrected_track = [];
for i = 1:length(all_t)
    
    %Select all entries with time all_t(i)
    IDs = find(cat_track(:,1) == all_t(i));
    
    %If more than one...
    if(length(IDs)>1)
        %Get all intensity values of these points. 
        fit_ids = cat_track(IDs,2);
        int = cat(1,fits_all( fit_ids ).sum_int);
        [~,max_id] = max( int );
        keep_entry = IDs(max_id);
        corrected_track = [corrected_track; cat_track(keep_entry,:)];
    else
        corrected_track = [corrected_track; cat_track(IDs,:)];
    end
end
        
final_fits = fits_all( corrected_track(:,2) );
corrected_track(:,2) = [1:size(corrected_track,1)]';

%Update obj info.
obj.spot_tracks{1} = corrected_track;
obj.results = final_fits;

%% Now get an ROI polygon for the background. 
%Imshow3D. 
figure(2);
imshow3D(IMG);

%Ask for ROI. 
h = impoly('Closed',1);

%Get positions of polygon. 
vec = h.getPosition;

%Convert to mask. 
MASK = poly2mask(vec(:,1),vec(:,2),size(IMG,1),size(IMG,2));

%Get pixel idx. 
PixelIdx = find(MASK);

%Now make a FRAP object. 
frap = frap_analyzer( obj.exp_info );
roi.PixelIdx = PixelIdx;
roi.type     = 'polygon';
roi.img_type = 'img';

%Calculate mean trace with open image. 
for i = 1:T
    this_img = IMG(:,:,i);
    int_trace(i) = mean( this_img( PixelIdx ) );
end
roi.trace = int_trace;
frap = frap.addROI( roi );

%Now add frap object to obj. 
obj.frap = frap;
disp('finished.')
close(gcf);

end

function empty_structure = gen_fit_struct( N )

%Fields
empty_structure(N) = struct('pos',[],'mean_int',[],...
    'real_bg',[],'sum_int',[],'cell_id',[],'snr',[]);

end

%Get fits based on filtering detection
function fit = get_fits(stack, params)

%Decide whether to do 2D or 3D. 
switch params.fit_type
    
    case '2D'

        mip = max(stack,[],3);
        h = fspecial('log',[15,15],8); %,[5,5],params.log_sigma(1));
        img_filter = -imfilter(mip,h,'symmetric', 'conv');

        img_bw = img_filter >= params.thresh;
        stats = regionprops(img_bw,'Centroid','Area');
        A = num2cell( ones(length(stats),1));
        [stats.assignment] = A{:};    
        if(isempty(stats))
            fit = [];
        else
            fit = fit_this_frame_centroid_2D( mip, stats, params );
        end
        
    case '3D'
        
        img_filter = -log_filter_3D(stack,params.sigma,params.filter_size);
        img_bw = img_filter >= params.thresh;
        stats = regionprops(img_bw,'Centroid','Area');
        A = num2cell( ones(length(stats),1));
        [stats.assignment] = A{:};
        if(isempty(stats))
            fit = [];
        else
            fit = fit_this_frame_centroid(stack,stats,params,step);
        end
end
        

end