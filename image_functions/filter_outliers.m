function BW = filter_outliers(I, BW, OutlierThreshold, MeanFilterNeighborhood )
        %Estimate background value. 
        bg_img = I(~BW);
        mean_bg_value = mean(bg_img(:));
        
        %Find outliers using a high threshold. 
        high = I > OutlierThreshold;
        
        %Dilate 
        se = strel('rectangle', [MeanFilterNeighborhood(1), MeanFilterNeighborhood(1)]);
        dil = imdilate(high, se);
        
        %Collect region info. 
        stats = regionprops(logical(dil),'PixelIdxList','Image');
        
        %Change high-intensity regions to NaN or mean bg value? 
        all_idx = cat(1,stats.PixelIdxList);
        I(all_idx) = mean_bg_value;
        
        BW(all_idx) = 0;
end