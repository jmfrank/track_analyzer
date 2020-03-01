%% Watershed
function [BW, stats] = watershed_objects(BW, nconn_BW, AbsMaxVol, h_min_depth, h_min_conn, z_effect, px_size)

    %Connect regions, find objects that are are big enough to be doubles. 
    Con=bwconncomp(logical(BW), nconn_BW);
    stats = regionprops(Con,'Centroid','Area','PixelList','PixelIdxList');

    %New volume list
    volumes = [stats.Area]';   
    liar_list = find( volumes > AbsMaxVol);
    disp('Starting watershed')
    
    
    img_dim= length(size(BW));
    
    switch img_dim
        
        
        case 2

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

        case 3
            

            %Distance transform scales. 
            %scales = [px_size(1),px_size(2), px_size(3)*z_effect];
            scales = [1,1,z_effect];
            %Counter of new objects found by watershedding
            new_objects = 0;

            %Empty fused object imglabel
            fused_imgLabel = uint16(zeros(size(BW))); 

            %Loop over liars
            for i = 1:length(liar_list)

                %This liar idx
                idx = liar_list(i);

                %Sub_img containing this liar_volume. 
                Px = stats( idx ).PixelList;

                %Pixel ranges. 
                min_range = min(Px,[],1);
                max_range = max(Px,[],1);

                x_range = min_range(1):max_range(1);
                y_range = min_range(2):max_range(2);
                z_range = min_range(3):max_range(3);

                %Create a sub image of BW. 
                sub_bw = zeros(length(y_range),length(x_range),length(z_range));

                %Add in the ones from this blob. 
                X = Px(:,1) - x_range(1) + 1;
                Y = Px(:,2) - y_range(1) + 1;
                Z = Px(:,3) - z_range(1) + 1;
                ind = sub2ind(size(sub_bw),Y,X,Z);
                sub_bw(ind) = 1;

                %Distance transform (faster without using gpu when image is
                %small. 
                imgDist = -bwdistsc(~sub_bw,scales);

                %Smoothing. Med filter is actually important here! 
                imgDist = medfilt3(imgDist,[5,5,1]);    

                %Get seeds  %%ORIGINAL params were 0.7,6. With medfilt3 = 5,5,5. 
                mask = imextendedmin(imgDist,h_min_depth,h_min_conn); %Seems like smaller neighborhood works better?
                imgDist = imimposemin(imgDist,mask);
                imgDist(~sub_bw) = -inf;

                %Perform marked watershed
                imgLabel = watershed(imgDist);
                %Replace background again?
                imgLabel(~sub_bw) = NaN;

                %Putting objects into new compiled BW. Figure out how many
                %objects exist after watershedding. Infer from the values of
                %imgLabel. imglabel=0 or 1 is background or envelope. 

                %Number of objects
                num_objects = double(max(imgLabel(:))-1 ); %Subtract 1 cause 0 and 1 are nothing. 

                %Ignore non-object pixels. 
                sel = imgLabel > 1;
                imgLabel(~sel) = 0;

                %Shift values by the on going counter. 
                imgLabel(sel) = imgLabel(sel) + new_objects;

                %Now place the imgLabel into the fused BW img back where it
                %came from. 
                fused_imgLabel(y_range, x_range, z_range ) = fused_imgLabel(y_range, x_range, z_range ) + uint16(imgLabel);

                %Adjust counter. 
                new_objects = new_objects + num_objects;
            end

            %Keep non liar stats
            keep = setdiff([1:length(stats)],liar_list);
            stats = stats(keep); 

            %Collect stats on watershedded regions. 
            stats_fused = regionprops(fused_imgLabel,'Centroid','Area','PixelList','PixelIdxList');
            stats_fused = stats_fused(2:end); %Ignore first entry. 
            %Append to stats. 
            stats = [stats; stats_fused];
            %Now we need to adjust BW... 
            BW = zeros(size(BW));
            %Add in new regions
            new_idx = cat(1,stats.PixelIdxList);
            BW(new_idx) = 1;
            
            disp('Finished watershed')

    end
end