%Go through frame_obj data and make tracks out of the nascent spot data. 

function obj = mk_spot_tracks(obj)


debug = 0;

%Clear out the any old data
obj.spot_tracks = cell(1,length(obj.tracks));


Z = obj.exp_info.z_frames;
T = obj.exp_info.t_frames;

%Get segmentation files 
    seg_files = dir([obj.exp_info.nuc_seg_dir,'*.mat']);

%Preload segfiles for faster access. 
    for f = 1:length(seg_files)
        %File name
        f_name = [obj.exp_info.nuc_seg_dir,seg_files(f).name];
        data(f) = load(f_name,'frame_obj'); 
    end
    
%The localization data is stored in data(i).frame_obj.results. 
%The variables are: {'Cellid','sigmaXY','sigmaZ','Y','X','Z','Int','BG','D2P'}

%Now we need to loop over the cell tracks, then make molecule tracks from this. 
%12-29-17: this script assumes that there's only 1 relevant spot per cell per time.
%We already segmented spots assuming the 1 spot rule (taking max
%intensity). 
for i = 1 : length( obj.tracks )
    %Get time vector
    time_vec = obj.tracks{i}(:,1);
    
    %First collect the cell id vector from this track (column 2). 
    cell_ids = obj.tracks{i}(:,2); 

    %Now go through time, and determine if cell of interest has a nascent
    %spot
    track = []; %<empty track, build it up
    for j = 1:length( time_vec )
        %current frame
        t = time_vec(j);
        %current cell id
        c = cell_ids(j);
        %Current frame_obj
        frame_obj = data( t ).frame_obj;
        %Check if there are spot tracking results for this frame
        if( ~isfield(frame_obj,'results') )
            continue
        end
        
        %Find if cell id 'c' is in frame_obj.results 
        sel = find( frame_obj.results(:,1) == c );
        
        %If sel is empty, this time point is empty
        if( isempty( sel ) )
            result = NaN(1,8);
            
            track = [track; c, t, 0, result];
        %If not, then append track data. Columns are: 
        % {cell id, frame, yes or no spot, 'sigmaXY','sigmaZ','Y','X','Z','Int','BG','D2P'}
        elseif( length(sel) == 1)
            %Grad result data
            result = frame_obj.results(sel,:);
            %Append to track
            track = [track; c, t, 1, result(2:end)];
        end

        
    end
    
    %Now add this track to the obj
    obj.spot_tracks{i} = track; 
end



end


