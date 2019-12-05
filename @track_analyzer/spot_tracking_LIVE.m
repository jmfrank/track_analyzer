 %Go through the raw LSM file and perform the spot tracking in 2D or 3D. 
 %This version uses a LoG filter (same as data-doctor routine) to detect
 %spots. This should work much better than DoG. 

 
 %6-7-18: changing to deal with 3D segmentation.
 %6-29-18: developed a way to locally threshold cells based using
 %percentile of intensities for each nucleus. 
 
function obj = spot_tracking_LIVE(obj, params, step, frames)

step = default_step( step )

% Image channel
channel_str = ['seg_channel_',pad(num2str(params.seg_channel),2,'left','0')];

%Generate reader. FOR NOW, assume we are looking in series 1. 
[reader,X,Y,Z,C,T] = bfGetReader(obj.exp_info.img_file);
series = 1;

%Figure out if call to function included specific frames. 
if nargin < 4
    frames = 1:T;
else
    frames = frames;
end

%Get segmentation files 
seg_files = obj.get_frame_files();

% Flags. 
if isfield(obj.exp_info,'flagged')
    disp('Found flagged data')
    if isempty(obj.exp_info.flagged)
        flags = [];
    else
        flags = obj.exp_info.flagged.cells;
    end
else
    flags = [];
end

%LSM time series goes Z first, then time. For each time point, load the
%Z-stack at time t, then max project and segment. 
for t = frames
    %Get the images for this z-stack according to bioformats reader
    img = zeros(Y,X,Z);
    for i = 1:Z
        this_plane = reader.getIndex(i-1,params.seg_channel-1,t-1)+1;
        img(:,:,i) = bfGetPlane(reader,this_plane);
    end
    
    %Load frame obj to get segmented cells. 
    load(seg_files{t}, 'frame_obj');
    
    % Deal with flags.
    indices = 1:length(frame_obj.(channel_str).centroids);
    if ~isempty(flags)
        ignore = flags(flags(:,1)==t,2);
        indices = setdiff(indices,ignore);
    end
        
    %Optional smoothing step.
    if step.smooth_img
       
        img = imgaussfilt3(img, params.smoothing_sigma);
    end
    
    %% Subtract local background? 
    if step.subtract_local_background 
        
        % Using the nuclear masks, subtract the median/mean intensity value
        % from each cell nuclei. This will prevent large gradients at
        % nuclear edge. 
        img = subtract_local_backgroud(img, frame_obj.(channel_str), indices) ;
    end
    
    if step.filter=='3D'
        
        %New LoG filter in 3D. Using fast implementation.
        img_filter = -log_filter_3D(img, params.log_sigma);
        
    elseif step.filter == '2D'
        
        img_filter = -log_filter_2D(img, params.log_sigma);
        
    end
    

    %Adding a local thresholding algorithm. Stats contains centroids and cell assignments.  
    stats = local_threshold_cells( img, img_filter, frame_obj, indices, params);
        
    %Check if there were any ROI's
    if( length( stats )==0 )
        %continue
    end
    
    % Merge duplicates (i.e. spots that are too close. 
    if step.merge_duplicates
        
       stats = merge_duplicates(stats,params);
                
    end
    
    % Filter out small or large regions. fit_this_frame_centroid
    if step.sizing
        areas = [stats.Area];
        sel = areas >= params.sizes(1) &  areas <= params.sizes(2);
        stats = stats(sel);
    end
    
    if step.debug
        figure(10)
        clf
        imshow3D(img_filter);
        hold on
        %Concat xy coordinates of centroids
        centroids = cat(1,stats.Centroid);
        centroids = centroids(:,[1,2]);
        for i = 1:length(stats)
            h = viscircles( centroids(i,:), 15);
            h.Children(1).UserData = i;           
        end
        
        cell_centroids = cat(1,frame_obj.(channel_str).centroids);
        for i = 1:length(cell_centroids)
            %text(cell_centroids{i}(1),cell_centroids{i}(2),num2str(i),'fontsize',20,'Color','w')
        end
        
        colormap gray
    end
    
    %Decide which fitting processes to perform
    switch step.Fit
        
        case '3D'
           
            %Perform fitting in 3D
            fit = fit_this_frame( img, stats, params, step );

            
        case '2D'
    
            %Use maximum projection image or guassial filtered image
            fit = fit_this_frame_mip( mip, stats, frame_obj, params, step );
            
        case 'Centroid'
            
            fit = fit_this_frame_centroid( img, stats, params);
        
        case 'Arbitrary'
            
            fit = fit_this_frame_arbitrary( img, stats, frame_obj.(channel_str).BW, params);
            
    end
    
    %% Remove fits if necessary.
    if step.max_fits_per_nucleus
        new_fits = [];
        cell_ids = [fit.cell_id];
        
        for i = unique(cell_ids)
            these_fits = cell_ids == i;
            ints = [fits(these_fits).sum_int];

            % Find by top intensity?
            I = [fit.sum_int];
            [Is,idx] = sort(I,'descend');
            
            if length(Is) > params.max_fits
                new_fits=[new_fits,fits(idx(1:params.max_fits))];
            end
            
        end
    end
    
    
    %Append fitting results to frame_obj
    frame_obj.(channel_str).fit = fit;
    save(seg_files{t}, 'frame_obj','-append');
    
    %struct2table(fit) %,'VariableNames',{'Cellid','sigmaXY','sigmaZ','Y','X','Z','Int','BG','D2P'})
    display(['Completed frame: ',num2str(t)])
end

reader.close();

end


%% 

function img = subtract_local_backgroud(img, frame_obj, indices)
    
    % Loop over cell nuclei in BW.     
    for i = 1:length(indices)
        
        bw = false(size(frame_obj.BW));
        
        bw( frame_obj.PixelIdxList{indices(i)} ) = 1;
        
        mean_val = mean( img(bw));
        %median_val = median( img(bw));
        
        img(bw) = img(bw)-mean_val;
        
    end
    
    % Set the background to zero?
    bg = ~frame_obj.BW;
    img(bg) = 0;
    
end
    

%% Merging centroids that are very close into one defined region. 
function stats = merge_duplicates(stats, params)


    D = pdist2( cat(1,stats.Centroid), cat(1,stats.Centroid) );

    these_too_close = triu(D) < params.merge_distance & triu(D) > 0;

    %Loop over and merge. 
    old_stats = stats;
    stats(length(old_stats))=struct('Centroid',[],'PixelIdxList',[],'assignment',[],'Area',[]);
    c=0;
    bad_list = [];
    for i = 1:size(these_too_close,1)
        [idx] = find(these_too_close(i,:));

        if any(i==bad_list)
            continue
        end

        if ~isempty(idx)

            % Iteratively check if each idx also have neighbors. 
            all_idx = idx;
            curr_idx = idx;
            pass = 0;
            while pass == 0
                sum_more = [];

                for j = 1:length(curr_idx)
                    sum_more = [sum_more,find(these_too_close(curr_idx(j),:))];
                end

                if isempty(sum_more)
                    pass = 1;
                else
                    all_idx = union(all_idx,sum_more);
                    curr_idx = sum_more;
                end
            end
            idx = all_idx;
            c=c+1;

            % Take mean position of all centroids. 
            stats(c).Centroid =  mean(cat(1,old_stats([i,idx]).Centroid));
            stats(c).PixelIdxList = unique(cat(1,old_stats([i,idx]).PixelIdxList));
            stats(c).assignment = old_stats(i).assignment;
            stats(c).Area = length(stats(i).PixelIdxList);
            bad_list = [bad_list, idx];
        else
            c=c+1;
            stats(c) = old_stats(i);
        end
    end
    %Remove empties. 
    stats = stats(1:c); 
end

%% defaults. 
function step = default_step( step )

%List of all default parameters. 
dstep.max_fits_per_nucleus=0;
dstep.debug=0;

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