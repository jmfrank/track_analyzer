%% Watershed
function [BW, stats] = watershed_objects(BW, nconn_BW, AbsMaxVol, h_min_depth, h_min_conn)

    %Connect regions, find objects that are are big enough to be doubles. 
    Con=bwconncomp(logical(BW), nconn_BW);
    stats = regionprops(Con,'Centroid','Area','PixelList','PixelIdxList');
    numobj=numel(stats); %total number of segmented objects

    %New volume list
    volumes = [stats.Area]';   
    liar_list = find( volumes > AbsMaxVol);
    disp('Starting watershed')
    %Make bw of clusters only
    bw_clust = zeros(size(BW));
    idx = cat(1,stats(liar_list).PixelIdxList);
    bw_clust(idx) = 1;
    
    %distance transform
    imgDist = -bwdist(~bw_clust);
    
    %Smoothing
    imgDist = medfilt2(imgDist,[3,3]);
    
    %Get seeds 
    mask = imextendedmin(imgDist, h_min_depth, h_min_conn); 
    
    %Seems like smaller neighborhood works better?
    imgDist = imimposemin(imgDist,mask);
    imgDist(~bw_clust) = -Inf;
    
    %Perform marked watershed
    imgLabel = watershed(imgDist);
    
    %Remove background again.
    imgLabel(~bw_clust) = NaN;
    stats_fused = regionprops(imgLabel,'Centroid','Area','PixelList','PixelIdxList');
    
    %Skip first entry. 
    stats_fused = stats_fused(2:end);

    %Now get rid of original regions that were fused. 
    keep = setdiff([1:length(stats)],liar_list);
    stats = stats(keep);
    
    %Now append the fused nuclei stats. 
    stats = [stats;stats_fused];
    
    %Now we need to adjust BW... 
    BW = zeros(size(BW));
    
    %Add in new regions
    new_idx = cat(1,stats.PixelIdxList);
    BW(new_idx) = 1;
    disp('Finished watershed')
end