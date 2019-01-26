%Combine data from multiple image series.

function combined = combine_experiments(info) %csv_file, exp_id, info_list, in_loc)

file_list = {};
img_list = {};
maxp_list = {};
time_series = [];
delT = 0;

%Loop over exp_ids
for i = 1:length(info.exp_id)
    info.exp_id(i)
    this_info.csv_file = info.csv_file;
    this_info.exp_id = info.exp_id(i);
    
    %Check location. 
    if isfield(info,'loc'); this_info.loc=info.loc; end;
    
    obj=get_exp( this_info );
    
    seg_files = obj.get_frame_files;
    
    time_series = [time_series, delT + obj.exp_info.time_series(1:length(seg_files))];
    if length(time_series) == 1
        %Ask user for input. 
        disp('Single time-point image:')
        disp(obj.exp_info.img_file);
        prompt = 'Enter the time till next experiment (in seconds):';
        delT = input(prompt);

    else
        delT = time_series(end) - time_series(end-1) + time_series(end);
    end
    
    file_list = [file_list; seg_files];
    
    img_list{i} = obj.exp_info.img_file;
    maxp_list{i} = obj.exp_info.max_p_img;
    
end

%Building up new info for combined data. 
new_info = obj.exp_info;
new_info = rmfield(new_info,{'sub_dir','img_file','track_file','nuc_seg_dir'});
new_info.time_series = time_series;
new_info.seg_files = file_list;
new_info.combined_experiments = [info.exp_id];
new_info.img_file = img_list;
new_info.max_p_img = maxp_list;

%Make new track_analyzer. 
this_info.csv_file=info.csv_file;
this_info.exp_id = info.out_id;
exp_info = get_exp_info(this_info);
combined = track_analyzer( exp_info );
combined = combined.update_exp_info( new_info );
combined.save;

%By default, start tracking, adding data
try
    params.max_dist= obj.exp_info.max_dist;
    combined = combined.track_cells(params);
    combined.save;
    disp('tracked successfully');
catch ME
    warning('could not track data');
    rethrow(ME)
end    
    
try
    combined = combined.add_nuc_cyto_signal_calcs;
    disp('added nuc/cyto data');
    combined.save;

catch ME
    warning('could not add N/C data');
    rethrow(ME)
end

end



