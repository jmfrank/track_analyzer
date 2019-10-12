%Tracking segmented objects based on centroids. 
%6-17-18: Updating to make sure it's compatible with 3D tracking. 
%9-4-18: Need to fix bug where there's no cells in the first frame. 
%11-10-18: Function looks for flagged data to ignore flagged cells. 
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

%% Look for flagged cells. 
if isfield( obj.exp_info,'flagged')
    
    flags = obj.exp_info.flagged.cells;
    if isempty(flags)
        flags=[0,0];
    end
else
    
    
    flags = [0,0]; %<<< dummy variable
    
end


%% Initialize tracks. Loop through frame_obj's until we see the first
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
    old_ids = [1:length(data.frame_obj.centroids)];
    
    % Inialize old_centroids
    old_centroids = cell2mat(data.frame_obj.centroids(old_ids)');
    old_centroid_loc = zeros(1,length(data.frame_obj.centroids));

    % Remove flagged data. 
    old_ids = setdiff(old_ids, flags( flags(:,1) == start_frame, 2));
    for i = 1:length(old_ids)
       obj.tracks{i}   = [start_frame, old_ids(i), old_centroids(i,:)];
       old_centroid_loc(old_ids(i)) = i;
       track_counter=track_counter+1;
    end
    

    %Now continue for all time points after start_frame. 
    disp_str = '';
    for i = start_frame+1:length(seg_files)
        
        fprintf( repmat('\b',[1,length(disp_str)+1]))
        disp_str = ['Analyzing frame: ',num2str(i)];
        disp(disp_str)
        
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
            %Ignore flagged data. 
            new_ids = [1:size(new_centroids,1)];
            new_ids = setdiff(new_ids,flags(flags(:,1)==i,2));
            new_centroids = new_centroids(new_ids,:);
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

                % ASSIGNMENT. output of munkres gives the assignment vector. For each column (
                % first frame centroids), it gives the optimal index for the
                % new_centroids
                [assignment, cost] = munkres(costMat);

                %The index returned corresponds to the new_centroid 
                %First go over non-zero assignments. 
                a_idx = find(assignment > 0);

                %Clear new_centroid_loc before going over tracks. allocate as zeros
                %- don't know the actual size....
                new_centroid_loc = [];
        
                %Loop over these indices
                for j =1:length(a_idx)

                    %Get centroid idx
                    new_centroid_id = new_ids(assignment(a_idx(j)));
                    new_centroid = new_centroids(assignment(a_idx(j)),:);
                    
                    %Get the track id
                    track_id         = old_centroid_loc( old_ids(a_idx(j)) );
                    
                    %Append the track with next centroid id/location
                    obj.tracks{track_id}   = [obj.tracks{track_id}; i, new_centroid_id, new_centroid];

                    %Need to update new_centroid_loc which keeps track of where
                    %centroids are placed. 
                    new_centroid_loc(new_centroid_id) = track_id;

                end
            %These are the centroids NOT assigned
            row_idx_new = setdiff(1:size(new_centroids,1), assignment(a_idx));
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
        if ~isempty(row_idx_new)
            
            for j = 1:length(row_idx_new)
                track_counter=track_counter+1;
                %Append new tracks at the end of the obj.tracks
                obj.tracks{track_counter}    = [i, new_ids(row_idx_new(j)), new_centroids(row_idx_new(j),:)];
                %Need to update new_centroid_loc which keeps track of where
                %centroids are placed. 
                new_centroid_loc(new_ids(row_idx_new(j))) = track_counter;           
            end
            
        end
        
        %Need to keep track of the track id location of every centroid to
        %match on next iteration! 
        old_centroid_loc = new_centroid_loc;
        %Now re-define the old_centroids 
        old_centroids=new_centroids;
        old_ids = new_ids;
        
    end
    
    %Remove empty tracks. 
    obj.tracks=obj.tracks(1:track_counter);
    
    %Remove references to flagged tracks. Flagged cells is still relevant? 
    try
        obj.exp_info.flagged = rmfield(obj.exp_info.flagged,'tracks');
    end
    
    if(debug);set(gca,'Ydir','reverse');end
end
