%Go through frame_obj data and make tracks out of the nascent spot data. 

function obj = make_spot_tracks(obj, seg_files)

debug = 0;

%Clear out the any old data
obj.spot_tracks = cell(1,length(obj.tracks));

%Get segmentation files 
    if(nargin==1)
        seg_files = obj.get_frame_files;
    end
    
    %Preload segfiles for faster access. 
    for f = 1:length(seg_files)
        %File name
        data(f) = load(seg_files{f},'frame_obj'); 
    end
    
    %Now append all fit data to results structure
    result = [];
    frames = [];
    for f = 1:length(data)
        these_fits = data(f).frame_obj.fit;
        %Check if there were any fits
        tmp_int = cat(1,these_fits.sum_int);
        if( ~isempty( tmp_int ) )
            result = [result, these_fits];
            frames = [frames; f.*ones(length(these_fits),1)];
        end
    end
       
    %Add result to object
    obj.results = result;
    
    %We'll also need the vector of cell ids for all localizations
    cell_ids_all = cat(1,result.cell_id);

%Now we need to loop over the cell tracks, then make molecule tracks from this. 
%12-29-17: this script assumes that there's only 1 relevant spot per cell per time.
pos_time_points = [];

for i = 1 : length( obj.tracks )
    %Get time vector (actually, this is the frame number)
    time_vec = obj.tracks{i}(:,1);
    
    %First collect the cell id vector from this track (column 2). 
    cell_ids = obj.tracks{i}(:,2); 

    %Now go through time, and determine if cell of interest has a nascent
    %spot
    count = 0;
    %Make a track variable. Will be column vectors of frame value and
    %localization id. 
    track = []; 
    
    %Loop over times
    for j = 1:length( time_vec )
        %current frame
        t = time_vec(j);
        %current cell id
        c = cell_ids(j);
                     
        %Get the index when frames is t and cell_id == c. 
        these_localization_ids = find(frames == t & cell_ids_all  == c);
        
        %Check if one or more localizations at this time. 
        if( length(these_localization_ids) == 1)
            %Add to counter
            count = count + 1;
            %Append track variable with frame value and id
            track(count,:) = [t,these_localization_ids];
            
        elseif( length(these_localization_ids) >1)
            
            %Find brightest spot. 
            count = count + 1;
            these_results = result(these_localization_ids);
            [~,max_id] = max( cat(1,these_results.sum_int));
            track(count,:) = [t,these_localization_ids( max_id )];
        end
    end
    
    %Now add this track to the obj
    obj.spot_tracks{i} = track;     
    display(['Finished track: ',num2str(i)])
end



end


