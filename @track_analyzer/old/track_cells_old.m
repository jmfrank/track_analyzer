%Function tracks segmented cells

function obj = track_cells(obj, params)

%Get main directory
main_dir = obj.exp_info.sub_dir; 
img_dir  = obj.exp_info.nuc_seg_dir
files = dir([img_dir,'*.mat'])
[img_dir,files(1).name];
track_file = [main_dir,'track_data.mat'];

debug = 0;

%Initialize tracks
data{1} = load([img_dir,files(1).name],'frame_obj');

%Assign tracks
    %The track object is a vector of indices. The value of the vector at a certain position corresponds to the index
    %of the contour for the time point corresponding to that position.
    %tracks have the format: [time point, centroid id from that time point,
    %centroid position x, y]. 
    
    ids = [1:length(data{1}.frame_obj.centroids)];
    for i = 1:length(ids)
       obj.tracks{i}   = [1,ids(i),data{1}.frame_obj.centroids{i}];
    end
    %Inialize old_centroids
    old_centroids = cell2mat(data{1}.frame_obj.centroids(:));
    old_centroid_loc = get_centroid_location(obj);
    %Add img files to data
    [a,b,c] = fileparts([obj.exp_info.nuc_seg_dir,files(1).name]);    
    obj.img_files{1}  = [a,'/',b,'.tif'];
    %Add cells to data? 
    obj.cells{1}  = data{1}.frame_obj.refined_cells;
    
    %Now continue for all time points
    for i = 2:length(files)
        display(['Analyzing frame: ',num2str(i)])
        
        %Load frame data
        data{i} = load([obj.exp_info.nuc_seg_dir,files(i).name],'frame_obj');
        %Add image files
        [a,b,c] = fileparts([obj.exp_info.nuc_seg_dir,files(i).name]);
        obj.img_files{i}  = [a,'/',b,'.tif'];
        %Add cell contours
        obj.cells{i}  = data{i}.frame_obj.refined_cells;
        img_file = obj.img_files{i};
        
        %For each centroid part of frame i, see which centroid is closest
        new_centroids = cell2mat(data{i}.frame_obj.centroids(:)); 
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
        % first frame centroids), it gives the optimal column index for the
        % new_centroids

        [assignment,cost] = munkres(costMat);
        
        %Conveniently, the index returned corresponds to the track id.
        %First go over non-zero assignments first. 
        row_idx = find(assignment > 0);

        color_vec  = hsv(length(row_idx));
        color_vec  = color_vec(randperm(size(color_vec,1)),:);
        
        %Clear new_centroid_loc before going over tracks
        clear new_centroid_loc
        %Go over matched tracks
        for j =1:length(row_idx)
            centroid_loc = assignment(row_idx(j));
            track_id         = old_centroid_loc( row_idx(j) );
            
            obj.tracks{track_id}   = [obj.tracks{track_id}; i, centroid_loc, new_centroids(centroid_loc,:)];

            
            %Need to update new_centroid_loc which keeps track of where
            %centroids are placed. 
            new_centroid_loc(centroid_loc) = track_id;
            
            %tmp_data = obj.tracks{track_id}
            time_vec     = obj.tracks{track_id}(:,1);
            track_id_vec = obj.tracks{track_id}(:,2);
            
            curr_pos = find(time_vec == i);
            prev_pos = find(time_vec == i-1);
            
            curr_cell_id = track_id_vec(curr_pos);
            prev_cell_id = track_id_vec(prev_pos);
            
            %track_id_vec(i-1)
            %data{i-1}
            cell_1 = data{i-1}.frame_obj.refined_cells{prev_cell_id};
            cell_2 = data{i}.frame_obj.refined_cells{curr_cell_id};
            if(debug)
                if(track_id==117)
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
                
        %Now we need to start new tracks for those not assigned. Need
        %to find the numbers not included in row_idx! This only occurs
        %when the number of new centroids is less than old centroids,
        %unless we have a fancier assignment algorithm with a maximum
        %cost. FIX LATER. 

        if( sum(assignment > 0) < size(new_centroids,1) )
            %Find the indices not used in the assignment variable. These
            %are the un-assigned tracks, so make all of these new tracks! 
            row_idx_new = setdiff([1:size(new_centroids,1)],assignment(row_idx));
            %Start a counter for new tracks 
            n=length(obj.tracks)+1;
            for j = 1:length(row_idx_new)
                %Append new tracks at the end of the obj.tracks
                obj.tracks{n}    = [i, row_idx_new(j), new_centroids(row_idx_new(j),:)];
                %Need to update new_centroid_loc which keeps track of where
                %centroids are placed. 
                new_centroid_loc(row_idx_new(j)) = n;
            
                %Plot these cells as black
                %cell_1 = data{i}.frame_obj.refined_cells{    }
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
    
    set(gca,'Ydir','reverse');    
end


%Returns the track locations of all centroids for the last frame. 
function loc = get_centroid_location(obj)
    %tracks
    tracks = obj.tracks;
    for i = 1:length(tracks)
        id = tracks{i}(end,2);
        loc(id) = i;
    end
end
