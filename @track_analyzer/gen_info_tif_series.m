%Create a structure for saving relevant image parameters from LSM files

function obj = gen_info_tif_series(obj, n_channels, n_time_points, pixel_size)


INFO = imfinfo(obj.exp_info.img_file);

%Numer of Z <<need to make sure this works? 
obj.exp_info.z_frames = length(INFO)./ n_channels / n_time_points;
%Number of frames and vector of time points. Assume 1 frame at T=0 for
%now...
obj.exp_info.t_frames = n_time_points;
obj.exp_info.time_series = 0; %The vector of time stamps for the image series (seconds?)

%Number of channels (assume 1 for now...)
obj.exp_info.Channels = n_channels;


%Try reading some metadata
try
    %Generate reader
    reader = bfGetReader(obj.exp_info.img_file);
    omeMeta = reader.getMetadataStore();
    
    %See if pixel size is available
    pixel_size(1) = omeMeta.getPixelsPhysicalSizeX(0).value();
    pixel_size(2) = omeMeta.getPixelsPhysicalSizeY(0).value();
    pixel_size(3) = omeMeta.getPixelsPhysicalSizeZ(0).value();
    
    obj.exp_info.pixel_size = pixel_size;
    
catch
    
    %Just look for pixel_size variable as input
    obj.exp_info.pixel_size = pixel_size;
end




end




