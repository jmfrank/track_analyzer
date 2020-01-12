function [BW, stats] = filter_objects(BW, min_size, max_size)
    if nargin==2
        max_size=Inf;
    end
    
    %Collect stats. 
    stats = regionprops(BW,'Centroid','Area','PixelList','PixelIdxList');

    %Remove very low volumes. 
    volumes = [stats.Area]';
    sel = volumes >= min_size & volumes <= max_size;
    stats = stats(sel);

    %Rebuild BW
    BW = zeros(size(BW));
    idx = cat(1,stats.PixelIdxList);
    BW(idx) = 1;
    BW = logical(BW);
end