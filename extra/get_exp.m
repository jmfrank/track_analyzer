% Retrieve data obj, or generate one if it doesn't exist. 
% info structure must contain csv_file and exp_id (row in csv_file)

function track_obj = get_exp(info) 

%Retrieve experiment info. 
exp_info = get_exp_info(info);

if ~isfield( info, 'make_new')
    info.make_new=0;
end


% Add bioformats javapath. 
class_path = which("track_analyzer");
S = split(class_path,'/');
bf_path = fullfile('/',S{1:end-2},'image_functions','bfmatlab','bioformats_package.jar');
javapath = javaclasspath('-all');   % Get the Java classpath

if all(cellfun(@isempty, strfind(javapath, 'bioformats_package.jar')))  %#ok<*STRCLFH>;
    javaaddpath(bf_path, '-end');
    disp(['Adding "' bf_path '" to Matlab java path']);
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
    catch ME
        warning('Could not perform bio-formats reading on img_file');
        rethrow(ME)
    end
else
    %Load object. 
    load(exp_info.track_file);
end

%Run update function. 
track_obj = track_obj.update_exp_info(exp_info);

end
