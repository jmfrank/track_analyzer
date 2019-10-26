%Go through tracks and add ids that point to the nuc_cyto_data found in
%frame_obj.
function obj = add_nuc_cyto_signal_calcs(obj, params)

debug = 0;

%Clear out old data? 
obj.nuc_cyto_data = [];

%Get segmentation files 
seg_files = obj.get_frame_files;
    
%Figure out channel type (could be channel specific or old 'data' type. 
tmp = load(seg_files{1},'frame_obj');

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
        channels=cellstr('');
        c=0;
        for i = 1:3
            if isfield(tmp.frame_obj,['channel_',pad(num2str(i),2,'left','0')])
                c=c+1;
                channels{c} = ['channel_',pad(num2str(i),2,'left','0')];
            end
        end
    end
end

    data=struct();
    %Preload segfiles for faster access. 
    for f = 1:length(seg_files)
        tmp = load(seg_files{f},'frame_obj'); 
        for c=1:length(channels)
            try
                data(f).(channels{c}) = tmp.frame_obj.(channels{c});
            catch
                warning(['No cells found in frame: ',num2str(f)]);
            end
        end
    end
    
    
    %Initialize empty data structure. 
    for c = 1:length(channels)
        nuc_cyto_data.(channels{c}) = gen_data_struct(10000);
    end
    
    %Loop over frames. Check that each channel exists, append. 
    frames = [];
    count=0;
    for f = 1:length(data)
        
        %Loop over channels. 
        for c = 1:length(channels)
            
            if ~isfield(data(f),channels{c})
                these_data=[];
                continue
            end
            
            these_data = data(f).(channels{c});
            new_length=length(these_data);
            nuc_cyto_data.(channels{c})(count+1:count+new_length) = these_data;
        end
        
        %Counter should be outside of channel loop. 
        count=count+new_length;

        %Because the channel data is 1-for-1, just need one frames vector. 
        if ~isempty(these_data)
            frames = [frames; f.*ones(length(these_data),1)];
        end
    end
    
    %Remove empty entries. 
    for c = 1:length(channels)
        nuc_cyto_data.(channels{c}) = nuc_cyto_data.(channels{c})(1:count);
    end
    
    %Add result to object
    obj.nuc_cyto_data = nuc_cyto_data;
    
    %We'll also need the vector of cell ids for all localizations
    cell_ids_all = cat(1,nuc_cyto_data.(channels{1}).cell_id);
    
    %Now loop over tracks, adding the ID to nuc_cyto_data as column 5 
    for i = 1 : length( obj.tracks )
        %Get time vector (actually, this is the frame number)
        time_vec = obj.tracks{i}(:,1);

        %First collect the cell id vector from this track (column 2). 
        cell_ids = obj.tracks{i}(:,2); 
        
        count = 0;
        %Make a track variable for vector of n/c data ids. 
        track = NaN(length(time_vec),1); 

        %Loop over times
        for j = 1:length( time_vec )
            %current frame
            t = time_vec(j);
            %current cell id
            c = cell_ids(j);

            %Get the index when frames is t and cell_id == c. 
            track(j) = find(frames == t & cell_ids_all  == c);
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