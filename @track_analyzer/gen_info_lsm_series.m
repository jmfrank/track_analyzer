%Create a structure for saving relevant image parameters from LSM files

function obj = gen_info_lsm_series(obj)

%LSM file. Can be time series or single time point. 
lsm_file = obj.exp_info.img_file;

%LSM reader
info = lsminfo(lsm_file);

%Get pixel size
obj.exp_info.pixel_size = info.VOXELSIZES(1)*10^6; %um
%Get frame to frame interval
obj.exp_info.frame_time = info.TimeInterval; %average time interval in milliseconds
%Get vector of frame times
obj.exp_info.time_series = info.TimeStamps.TimeStamps; %The vector of time stamps for the image series (seconds?)
%Number of time frames in this file. 
obj.exp_info.t_frames = length(obj.exp_info.time_series);
%Numer of Z <<need to make sure this works? 
obj.exp_info.z_frames = info.DimensionZ;

%Number of channels
obj.exp_info.Channels = info.ScanInfo.IMAGES_NUMBER_CHANNELS;
end




