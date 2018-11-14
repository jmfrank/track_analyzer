%Read experiment information for LSM series. Used for Nascent rna locus
%tracking. 

function track_obj = get_exp(info) %csv_file, exp_id, info_list, in_loc)

%Get experiment info. 
exp_info = get_exp_info(info);

if ~isfield( info, 'make_new')
    info.make_new=0;
end

%Check if experiment already exists. 
if ~exist(exp_info.track_file,'file') | info.make_new
    track_obj=track_analyzer(exp_info);
    disp('Making new track analyzer.');
    
    %Try bioformats reading. 
    try
        if isfield(info,'pixel_size') & isfield(info,'time_interval')
            track_obj = track_obj.gen_info_bioformats(info.pixel_size,info.time_interval);
        else
             track_obj = track_obj.gen_info_bioformats();
        end
    catch
        warning('Could not perform bio-formats reading on img_file');
    end
else
    %Load object. 
    load(exp_info.track_file);
end

%Run update function. 
track_obj = track_obj.update_exp_info(exp_info);

end