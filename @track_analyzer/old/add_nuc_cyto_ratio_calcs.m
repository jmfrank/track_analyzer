%Go through tracks and add ids that point to the nuc_cyto_data found in
%frame_obj.


function obj = add_nuc_cyto_ratio_calcs(obj)

debug = 0;

%Clear out old data? 
obj.nuc_cyto_data = [];


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
    
    
    %Now append all fit data to results structure
    nuc_cyto_data = gen_data_struct( 1 );
    frames = [];
    for f = 1:length(data)
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
            end
        end

        %Now add this track to the obj
        obj.tracks{i}(:,5) = track;     
        display(['Finished track: ',num2str(i)])
    end

    
end



%% Generate an empty structure for keeping track of nuclear and cytoplasmic signals
function data = gen_data_struct( N )

%Fields
data(N) = struct('nuc_mean',[],'cyto_mean',[],'local_rho',[],'cell_id',[]);


end