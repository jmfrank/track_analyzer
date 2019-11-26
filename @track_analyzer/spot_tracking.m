%Perform spot detection in cell nuclei, 2D or 3D. Used for nascent transcription spots seen using FISH staining.
%Use spot_tracking_LIVE for live-cell imaging. 

function obj = spot_tracking(obj, params, step, frames)


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
    
    %% Subtract local background? 
    if step.subtract_local_background 
        
        % Using the nuclear masks, subtract the median/mean intensity value
        % from each cell nuclei. This will prevent large gradients at
        % nuclear edge. 
        img = subtract_local_backgroud(img, frame_obj.(channel_str)) ;
    end
    
    %New LoG filter in 3D. Using fast implementation.
    img_filter = -log_filter_3D(img, params.log_sigma);

    %Binarize image based on threshold. 
    im_threshold = img_filter >= params.threshold;


    %Now just multily by the cell nucleus mask
    %Check if BW is 2D or 3D. 
    if(numel(size(frame_obj.(channel_str).BW))==2)
        BW = repmat(frame_obj.BW,[1,1,size(im_threshold,3)]);
    else
        BW = frame_obj.(channel_str).BW;
    end
        
    im_threshold = im_threshold.*BW;
        
    %Collect stats on regions
    stats = regionprops(logical(im_threshold), 'Centroid', 'PixelIdxList', 'Area');
    
    %Assignment
    seg_dim = length(stats(1).Centroid);
    spot_centroids = cat(1,stats.Centroid);
    cell_centroids = cat(1,frame_obj.(channel_str).centroids{:});
    if seg_dim == 2 && length(size(BW))==2
        D = pdist2(spot_centroids, cell_centroids);
    elseif seg_dim ==2 && length(size(BW)) ==3
        cell_centroids = [cell_centroids, h/2*ones(size(cell_centroids,1),1)];
        D = pdist2(spot_centroids, cell_centroids);
    elseif seg_dim ==3 && length(size(BW))==3
        D = pdist2(spot_centroids, cell_centroids);
    end

    [~,assignment] = min(D,[],2);

    A = num2cell(assignment);
    [stats.assignment] = A{:};
    
    % Merge duplicates (i.e. spots that are too close. 
    if step.merge_duplicates
        
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


    %Check if there were any empty ROI's
    %if( length( stats )==0 )
    %    continue
    %end
    
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
            
            fit = fit_this_frame_arbitrary( img, stats, params);
    end
    
    %Append fitting results to frame_obj
    frame_obj.(channel_str).fit = fit;
    save(seg_files{t}, 'frame_obj','-append');
    
    struct2table(fit) %,'VariableNames',{'Cellid','sigmaXY','sigmaZ','Y','X','Z','Int','BG','D2P'})
    display(['Completed frame: ',num2str(t)])
end



end

%% 

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


