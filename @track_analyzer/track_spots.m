%Tracking segmented objects based on centroids. 
%6-17-18: Updating to make sure it's compatible with 3D tracking. 
%7-12-18: this version is for tracking spot localizations. 
function obj = track_spots(obj, params, seg_files)

%Parse input. Use input file list if provided. 
if(nargin == 2)
    %Get segmentation files
    seg_files = obj.get_frame_files;
end

debug = 0;
%Initialize tracks
load(seg_files{1},'frame_obj');
data.centroids = cat(1,frame_obj.fit.pos);
%Clear the tracks field. Start with 100 tracks.
obj.spot_tracks = cell(100,1);

%Assign tracks
    %The track object is a vector of indices. The value of the vector at a certain position corresponds to the index
    %of the contour for the time point corresponding to that position.
    %tracks have the format: 
    %[time point, centroid id from that time point, centroid position x, y, z]. 
    
    %Initialize tracking by adding all objects in first frame as new
    %tracks.
    ids = [1:length(data.centroids)];
    for i = 1:length(ids)
       obj.spot_tracks{i}   = [1,ids(i),data.centroids(i,:)];
    end
    %Inialize old_centroids
    old_centroids = cell2mat(data.frame_obj.centroids(:));
    old_centroid_loc = get_centroid_location(obj);
    %Add img files to data
    [a,b,c] = fileparts(seg_files{1});    
    %Add cells to data
    try obj.cells{1}  = data.frame_obj.refined_cells; end
    
    %Now continue for all time points
    for i = 2:length(seg_files)
        display(['Analyzing frame: ',num2str(i)])
        
        %Load frame data
        data = load(seg_files{i},'frame_obj');
        %Add image files
        [a,b,c] = fileparts(seg_files{i});
        %Add cell contours
        try obj.cells{i}  = data.frame_obj.refined_cells;end
        
        %For each centroid part of frame i, see which centroid is closest
        new_centroids = cell2mat(data.frame_obj.centroids(:)); 
        %Create distance matrix
            % D = pdist2( X , Y ). D is a mx by my matrix with         
            %(i,j) = dist( x(i), y(j) )
            % In our case, x is the i-1 frame, and y is the i frame. 
            D = pdist2(old_centroids, new_centroids,'Euclidean');
            %Now set the pairs above the max_dist to Inf cost
            sel = D > params.max_dist;
            D(sel) = Inf;
            % For now, cost is just D^2. This exaggerates the cost of
            % distance
            costMat = D; %.^2;
        
        % ASSIGNMENT. output of munkres gives the assignment vector. For each row (
        % first frame centroids), it gives the optimal index for the
        % new_centroids
        [assignment,cost] = munkres(costMat);
        
        %Conveniently, the index returned corresponds to the track id.
        %First go over non-zero assignments. 
        row_idx = find(assignment > 0);

        %Clear new_centroid_loc before going over tracks. allocate as zeros
        %- don't know the actual size....
        clear new_centroid_loc
        
        %Loop over these indices
        for j =1:length(row_idx)
            
            %Get centroid idx
            centroid_loc = assignment(row_idx(j));
            %Get the track id
            track_id         = old_centroid_loc( row_idx(j) );
            %Append the track with next centroid id/location
            obj.spot_tracks{track_id}   = [obj.spot_tracks{track_id}; i, centroid_loc, new_centroids(centroid_loc,:)];

            %Need to update new_centroid_loc which keeps track of where
            %centroids are placed. 
            new_centroid_loc(centroid_loc) = track_id;
            

            %Debuging section
            if(debug)
                color_vec  = hsv(length(row_idx));
                color_vec  = color_vec(randperm(size(color_vec,1)),:);
        
                if(track_id==117)
                    time_vec     = obj.spot_tracks{track_id}(:,1);
                    track_id_vec = obj.spot_tracks{track_id}(:,2);

                    curr_pos = find(time_vec == i);
                    prev_pos = find(time_vec == i-1);


                    curr_cell_id = track_id_vec(curr_pos);
                    prev_cell_id = track_id_vec(prev_pos);
            
                    cell_1 = data{i-1}.frame_obj.refined_cells{prev_cell_id};
                    cell_2 = data.frame_obj.refined_cells{curr_cell_id};
            
                    c1 = plot(cell_1(:,1),cell_1(:,2),'color',color_vec(j,:),'linewidth',3);
                    hold on
                    c2 = plot(cell_2(:,1),cell_2(:,2),'color',color_vec(j,:),'linewidth',3);

                    %label(c1,['T:',num2str(i-1),' cID#',num2str(prev_cell_id),' tID#',num2str(track_id)])
                    %label(c2,['T:',num2str(i),' cID#',num2str(curr_cell_id),' tID#',num2str(track_id)])
                    hold off    
                    pause
                end
            end
        end
                
        %Now we need to start new tracks for centroids not assigned to previously existing tracks. Need
        %to find the numbers not included in row_idx! This only occurs
        %when the number of new centroids is less than old centroids,
        %unless we have a fancier assignment algorithm with a maximum
        %cost. 

        %Check if there are unassigned centroids
        if( sum(assignment > 0) < size(new_centroids,1) )
            %Find the indices not used in the assignment variable. These
            %are the un-assigned tracks, so make all of these new tracks! 
            row_idx_new = setdiff([1:size(new_centroids,1)],assignment(row_idx));
            %Start a counter for new tracks 
            n=length(obj.spot_tracks)+1;
            for j = 1:length(row_idx_new)
                %Append new tracks at the end of the obj.spot_tracks
                obj.spot_tracks{n}    = [i, row_idx_new(j), new_centroids(row_idx_new(j),:)];
                %Need to update new_centroid_loc which keeps track of where
                %centroids are placed. 
                new_centroid_loc(row_idx_new(j)) = n;
            
                %Plot these cells as black
                %cell_1 = data.frame_obj.refined_cells{    }
                %Increase count by one
                n=n+1;
            end
        end
        %pause(0.2)
        %Need to keep track of the track id location of every centroid to
        %match on next iteration! 
        old_centroid_loc = new_centroid_loc;
        %Now re-define the old_centroids 
        old_centroids=new_centroids;
    end
    
    if(debug);set(gca,'Ydir','reverse');end
end


%Returns the track locations of all centroids for the last frame. 
function loc = get_centroid_location(obj)
    %tracks
    tracks = obj.spot_tracks;
    for i = 1:length(tracks)
        id = tracks{i}(end,2);
        loc(id) = i;
    end
end
