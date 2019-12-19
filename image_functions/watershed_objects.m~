%% Merge 0bjects close together.
function stats = merge_objects(stats, merge_dist, AbsMaxVol)
        %Centroids. 
        all_ctrs = cat(1,stats.Centroid); 
        
        %Pairwise distance. 
        D = pdist2(all_ctrs,all_ctrs);
        
        %Select pair-wise distances less than threshold. 
        S = D < merge_dist;
        
        %Get pairs. 
        [a,b] = find( tril( S ) );
        
        %Entries to get rid of. 
        bad_list = false(1,length(stats));
        c=0;
        
        %Merge algorithm. Loop over pairs, merge if result isn't too large.
        for m = 1:length(a)
            id_a = a(m);
            id_b = b(m);
            w_a = length(stats(id_a).PixelIdxList);
            w_b = length(stats(id_b).PixelIdxList);
            merge_vol = w_a + w_b;
            
            if merg_vol < AbsMaxVol
                bad_list(id_a) = 1;
                bad_list(id_b) = 1;
                new_list = [stats(id_a).PixelIdxList; stats(id_b).PixelIdxList];
                new_centroid = w_a/merge_vol*stats(id_a).Centroid + w_b/merge_vol*stats(id_b).Centroid;
                new_pixel_list = [stats(id_a).PixelList; stats(id_b).PixelList];
                c=c+1;
                new_stats(c).Centroid=new_centroid;
                new_stats(c).PixelList = new_pixel_list;
                new_stats(c).PixelIdxList = new_list;
            end
        end
        
        %Now remove bad objects, and append new ones. 
        stats = stats(~bad_list);
        stats = [stats; new_stats];
end