%% Parse divisions from the groups structure generated by gui. 
function div = parse_divisions_nc( obj, params )


 %Check if there are divisions. 
if isfield( obj.exp_info, 'groups') 
    div = obj.exp_info.groups;
else
    disp('No divisions found')
    div=[];
    return
end

%Loop over divisions, analyze tracks / frames. 
for d = 1:length(div)


    %% Defining the parent and daughters. 
    
    %Check for reference frame. Should be closest time to
    %metaphase. 
    if isempty( div(d).marked_frame )
        error('missing division frame marker!');
    end

    %Frame at which division reference is marked. 
    div_frame = div(d).marked_frame; 
    
    
    %Div tracks. 
    track_ids = div(d).cell_tracks;

    %Collect frames for all tracks associated with this division. 
    frames = cell(length(track_ids) ,1);
    for t = 1:length(track_ids)
        frames{t} = obj.tracks{ track_ids(t) }(:,1);
    end

    %First frame of each cell track. 
    first_frame = cellfun(@(x) x(1), frames );

    %Check if all tracks occur on/after marked frame (no detectable parent). 
    idx_before = find( first_frame < div_frame );
    
    if ~isempty(idx_before)
        
        %Must be a parent cell. Should be only one index.   
        [~,min_idx] = min( first_frame( idx_before ) );
        
        div(d).first_parent_frame = obj.tracks{ track_ids( idx_before(min_idx))}(1,1);
        parent_track_idx = track_ids(idx_before);

    else
        
        %No parent cell available. 
        parent_track_idx = [];        

    end  
    
    %% Now we need to find all daughter tracks.
    
    %Loop over all tracks, find if any frames are after marked time point. 
    daughter_track_idx = [];
    for t = 1:length( track_ids )
        
        %All frames of this track
        these_frames = obj.tracks{ track_ids(t) }(:,1);
        
        if any( these_frames > div_frame )
            
            %Mark this track as a daughter. 
            daughter_track_idx = [daughter_track_idx, track_ids(t)];
        end
    end
    
    %If one or two daughters, easy. if more, then throw error - deal with
    %later?
    if length(daughter_track_idx) <= 2 
            
        %Good number of tracks. 
        daughter_tracks = obj.tracks( daughter_track_idx );
    else

        warning('more than 2 daughters!');
        daughter_tracks = obj.tracks( daughter_track_idx );

    end

    %Define the MIN of the last observable frame of each daughter. 
    last_frame_obserable = min( cellfun(@(x) x(end,1), daughter_tracks ));
    div(d).daughter_last_frame_observable = last_frame_obserable;

    
    %Add parent / daughter ids. 
    div(d).parent_track_idx = parent_track_idx;
    div(d).daughter_track_idx = daughter_track_idx;
    


    %% Make a division track. 
    %Extract marked time. 
    mark = div(d).marked_frame;

    %Get all frames for parent/daughter. 
    parent_frames = cellfun(@(x) x(:,1), obj.tracks( div(d).parent_track_idx ),'uniformoutput',0);
    daughter_frames = cellfun(@(x) x(:,1), obj.tracks( div(d).daughter_track_idx ),'uniformoutput',0);

    %Time vector wrt to cytokinesis. 
    frames = unique([cat(1,parent_frames{:});cat(1,daughter_frames{:})]) - mark;

    NC = zeros(size(frames));
    counts = zeros(size(frames));
    
    %Loop over parents, then daughters, take mean N/C over time. 
    %All cell_tracks involved. 
    IDS = unique([parent_track_idx,daughter_track_idx]);
    for i = 1:length(IDS)
        
        this_track = obj.tracks{IDS(i)};
        %Get signal data. 
        data = obj.get_track_data(IDS(i));

        %Ratio
        sig = (cat(1,data.nuc_mean))./ ( cat(1,data.cyto_mean));        
        
        %Intersection with 'frames'
        [~,idx_this, idx_all] = intersect(this_track(:,1)-mark,frames);
        
        %Add signal. 
        NC(idx_all) = NC(idx_all) + sig(idx_this);
        counts(idx_all) = counts(idx_all) + 1;
        
    end
    
    %Mean NC. 
    muNC = NC ./ counts;
    
    %Make d_track. 
    div(d).NC = [frames,muNC];
    
    
end


end