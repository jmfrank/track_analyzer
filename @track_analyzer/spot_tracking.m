%Perform spot detection in cell nuclei, 2D or 3D. Used for nascent transcription spots seen using FISH staining.
%Use spot_tracking_LIVE for live-cell imaging. 

function obj = spot_tracking(obj, params, step, frames)


% Image channel to segment
channel_str = ['channel_',pad(num2str(params.seg_channel),2,'left','0')];

% Image channel that contains cell masks. 
cell_channel_str = ['seg_channel_',pad(num2str(params.cell_channel),2,'left','0')];


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
seg_files = obj.get_frame_files;

%LSM time series goes Z first, then time. For each time point, load the
%Z-stack at time t, then max project and segment. 
for t = frames

    
    %Get the images for this z-stack according to bioformats reader
    img = zeros(Y,X,Z);
    for i = 1:Z
        this_plane = reader.getIndex(i-1,params.seg_channel-1,t-1)+1;
        img(:,:,i) = bfGetPlane(reader,this_plane);
    end

    %Load cell data and collect centroids. 
    load(seg_files{t}, 'frame_obj');
    
    %IMPORTANT: we often deal with stitched images. So we should fill in
    %zero-valued pixels. 
    %img = fill_edges(img);
    
    % intialize segmenter. 
    S = spot_segmenter(img, frame_obj.(cell_channel_str), params);
    seg_steps = obj.exp_info.steps.(channel_str).spots;
    for i = 1:length(seg_steps)
        
        disp(seg_steps(i).type);
        S.(seg_steps(i).type);
        
    end
    
    stats = S.stats;
    
    
    %Decide which fitting processes to perform
    switch step.Fit
        
        %3d gaussian. 
        case '3D'
           
            %Perform fitting in 3D
            fit = fit_this_frame( img, stats, params, step );
           
        %2d gaussian
        case '2D'
            
            %Max project
            mip = max(img,[],3);
            %Use maximum projection image or guassial filtered image
            fit = fit_this_frame_mip( mip, stats, frame_obj, params, step );
    
        case 'Centroid'
            
            %Simple centroid fitting scheme
            fit = fit_this_frame_centroid( img, stats, params);
        
        case 'Arbitrary'
            
            fit = fit_this_frame_arbitrary( img, stats, frame_obj.(cell_channel_str).BW, params);
    end
    
    %Append fitting results to frame_obj
    frame_obj.(['seg_',channel_str]).fit = fit;
    save(seg_files{t}, 'frame_obj','-append');
    
    struct2table(fit) %,'VariableNames',{'Cellid','sigmaXY','sigmaZ','Y','X','Z','Int','BG','D2P'})
    display(['Completed frame: ',num2str(t)])
end



end

%% Local background subtracking.
function img = subtract_local_backgroud(img, frame_obj)

    % Loop over cell nuclei in BW. 
    n = length(frame_obj.PixelIdxList);
    
    for i = 1:n
        
        bw = false(size(frame_obj.BW));
        
        bw( frame_obj.PixelIdxList{i} ) = 1;
        
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
