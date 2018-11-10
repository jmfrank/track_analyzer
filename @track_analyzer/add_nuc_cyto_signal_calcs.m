%Go through tracks and add ids that point to the nuc_cyto_data found in
%frame_obj.
function obj = add_nuc_cyto_signal_calcs(obj, params)

debug = 0;

%Clear out old data? 
obj.nuc_cyto_data = [];

%Get segmentation files 
seg_files = obj.get_frame_files;
    
%Figure out channel type (could be channel specific or old 'data' type. 
tmp = load(seg_files(f),'frame_obj');

if nargin==2
    channels = params.channels;
    %Check these channels exist. 
    for i = 1:length(channels)
        if ~isfield(tmp.frame_obj,channels{i})
            disp(['missing channel: ',channels{i}])
        end
    end     
else
    
    %Check field names. 
    if isfield(tmp.frame_obj,'data')
        channels{1}='data';
    else
        channels=cellstr('')
        for i = 1:3
            if isfield(tmp.frame_obj,['channel_',
                channels{




    %Preload segfiles for faster access. 
    for f = 1:length(seg_files)
        data(f) = load(seg_files{f},'frame_obj'); 
    end
    
    
    %Now append all fit data to results structure
    nuc_cyto_data = gen_data_struct( 1 );
    %Check if nuc-area exists << Get rid of this at some point. 
    if ~isfield( data(1).frame_obj.data,'nuc_area')
       nuc_cyto_data = rmfield(nuc_cyto_data,'nuc_area');
    end
    
    frames = [];
    for f = 1:length(data)
        if( ~isfield(data(f).frame_obj,'data'))
            continue
        end
        these_data = data(f).frame_obj.data;
        nuc_cyto_data = [nuc_cyto_data, these_data];
        frames = [frames; f.*ones(length(these_data),1)];
    end
    
    %Remove empty entry
    nuc_cyto_data = nuc_cyto_data(2:end);   
       
    %Add result to object
    obj.nuc_cyto_data = nuc_cyto_data;
    
    %We'll also need the vector of cell ids for all localizations
    cell_ids_all = cat(1,nuc_cyto_data.cell_id);

    
    %Now loop over tracks, adding the ID to nuc_cyto_data as column 5 
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
            this_localization_id = find(frames == t & cell_ids_all  == c);

            %Make sure sel is length 1
            if( length(this_localization_id) == 1)
                %Add to counter
                count = count + 1;
                %Append track variable with frame value and id
                track(count) = [this_localization_id];

            else
                %No localization found for this time in this cell
                disp('missing data?')
                [j,t,c]
            end
        end

        %Now add this track to the obj
        obj.tracks{i}(:,5) = track(:);     
    end

    
end



%% Generate an empty structure for keeping track of nuclear and cytoplasmic signals
function data = gen_data_struct( N )

%Fields
data(N) = struct('nuc_mean',[],'cyto_mean',[],'local_rho',[],'cell_id',[],'nuc_med',[],'cyto_med',[],'nuc_area',[]);


end