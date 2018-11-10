%Tracking segmented objects based on centroids. 
%6-17-18: Updating to make sure it's compatible with 3D tracking. 
%9-4-18: Need to fix bug where there's no cells in the first frame. 
function obj = track_cells(obj, varargin)

p = inputParser;
p.addRequired('params',@isstruct);
p.addParameter('seg_files','default');
p.addParameter('bad_frames',[]);
p.parse(varargin{:});

params = p.Results.params;
bad_frames = p.Results.bad_frames;
if strcmp(p.Results.seg_files,'default')
    %Get segmentation files
    seg_files = obj.get_frame_files;
else
    seg_files = p.Results.seg_files;
end

%Update exp_info with the 'max_dist' param. 
obj.exp_info.max_dist = params.max_dist;

debug = 0;

%Initialize tracks. Loop through frame_obj's until we see the first
%segmented cell. 
for i = 1:length(seg_files)
    
    if any(bad_frames==i)
        continue
    end
    
    D = load(seg_files{i},'frame_obj');
    %Check if there are cells. 
    if(isfield(D.frame_obj,'centroids') )
        if( ~isempty(D.frame_obj.centroids) )
            start_frame = i;
            data = load(seg_files{start_frame},'frame_obj');
            break
        end
    end
end

%Error if there were no cells in any frames.
if(~exist('start_frame','var'))
    error('No cells found in any frame. Stopping.');
end

%Clear the tracks field
obj.tracks = cell(1,5000);
track_counter=0;
%Assign tracks
    %The track object is a vector of indices. The value of the vector at a certain position corresponds to the index
    %of the contour for the time point corresponding to that position.
    %tracks have the format: 
    %[time point, centroid id from that time point, centroid position x, y, z]. 
    
    %Initialize tracking by adding all objects in start_frame frame as new
    %tracks.
    ids = [1:length(data.frame_obj.centroids)];
    for i = 1:length(ids)
       obj.tracks{i}   = [start_frame,ids(i),data.frame_obj.centroids{i}];
       track_counter=track_counter+1;
    end
    %Inialize old_centroids
    old_centroids = cell2mat(data.frame_obj.centroids(:));
    old_centroid_loc = get_centroid_location(obj,track_counter);

    %Now continue for all time points after start_frame. 
    for i = start_frame+1:length(seg_files)
        display(['Analyzing frame: ',num2str(i)])
        
        %Load frame data
        data = load(seg_files{i},'frame_obj');
        %Now check if there any centroids in this frame.
        if(~isfield(data.frame_obj,'centroids') || any(bad_frames==i) )
            %There's no new centroids. 
            empty_frame = 1;
            old_centroids = [];
            %Go to next frame. 
            continue
        else
            %For each centroid part of frame i, see which centroid is closest
            new_centroids = cell2mat(data.frame_obj.centroids(:)); 
            empty_frame = 0;
        end
        
        %If there were previous centroids, then proceed forward to connect
        %them. If not, then add all new tracks. 
        if(~isempty(old_centroids))

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
                    obj.tracks{track_id}   = [obj.tracks{track_id}; i, centroid_loc, new_centroids(centroid_loc,:)];

                    %Need to update new_centroid_loc which keeps track of where
                    %centroids are placed. 
                    new_centroid_loc(centroid_loc) = track_id;


                    %Debuging section
                    if(debug)
                        color_vec  = hsv(length(row_idx));
                        color_vec  = color_vec(randperm(size(color_vec,1)),:);

                        if(track_id==117)
                            time_vec     = obj.tracks{track_id}(:,1);
                            track_id_vec = obj.tracks{track_id}(:,2);

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
            %These are the centroids NOT assigned
            row_idx_new = setdiff([1:size(new_centroids,1)],assignment(row_idx));
        else
            %No old centroids, so all centroids are new. 
            row_idx_new = 1:size(new_centroids,1);

        end
                
        %Now we need to start new tracks for centroids not assigned to previously existing tracks. Need
        %to find the numbers not included in row_idx! This only occurs
        %when the number of new centroids is less than old centroids,
        %unless we have a fancier assignment algorithm with a maximum
        %cost. 

        %Check if there are unassigned centroids that will be new tracks. 
        if(~isempty(row_idx_new) )
            %Start a counter for new tracks 
            
            for j = 1:length(row_idx_new)
                track_counter=track_counter+1;
                %Append new tracks at the end of the obj.tracks
                obj.tracks{track_counter}    = [i, row_idx_new(j), new_centroids(row_idx_new(j),:)];
                %Need to update new_centroid_loc which keeps track of where
                %centroids are placed. 
                new_centroid_loc(row_idx_new(j)) = track_counter;           
            end
        end
        
        %Need to keep track of the track id location of every centroid to
        %match on next iteration! 
        old_centroid_loc = new_centroid_loc;
        %Now re-define the old_centroids 
        old_centroids=new_centroids;
    end
    
    %Remove empty tracks. 
    obj.tracks=obj.tracks(1:track_counter);
    
    %Remove references to flagged tracks. Flagged cells is still relevant? 
    try
        obj.exp_info = rmfield(obj.exp_info.flagged,'tracks');
    end
    
    if(debug);set(gca,'Ydir','reverse');end
end


%Returns the track locations of all centroids for the last frame. 
function loc = get_centroid_location(obj,c)
    %tracks
    tracks = obj.tracks;
    for i = 1:c
        id = tracks{i}(end,2);
        loc(id) = i;
    end
end
