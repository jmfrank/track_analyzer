%Create a structure for saving relevant image parameters from LSM files
function obj = gen_info_bioformats(obj, pixel_size, time_interval)


%Creater reader. Try 
reader = bfGetReader(obj.exp_info.img_file);
omeMeta = reader.getMetadataStore();


%Try to get physical size. Everything in microns...Might not be there. 
try
    x = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER);
    y = omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER);
    % check if z is empty
    if isempty(omeMeta.getPixelsPhysicalSizeZ(0))
        z = 0; 
    else
        z = omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER);
    end
    obj.exp_info.pixel_size = [double(x),double(y),double(z)];
catch
    
    %Need to  supply the pixel_size manually.
    obj.exp_info.pixel_size = pixel_size;
end

%Get number of z-planes
obj.exp_info.z_planes = reader.getSizeZ;
obj.exp_info.t_frames = reader.getSizeT;
obj.exp_info.channels = reader.getSizeC;
obj.exp_info.img_size = [reader.getSizeY,reader.getSizeX];

%Try to get frame time interval if possible
try 
    T = reader.getSizeT;
    if(T>1)

        %Find the first plane of second z-series. 
        this_plane = reader.getIndex(0,0,1);
        obj.exp_info.frame_time = double(omeMeta.getPlaneDeltaT(0, this_plane).value);
        %Loop over time points to get the time-series vector. All in seconds. 
        for i= 1:reader.getSizeT
            this_plane = reader.getIndex(0,0,i-1);
            time_vec(i) = double(omeMeta.getPlaneDeltaT(0,this_plane).value());
        end
    
        obj.exp_info.time_series = time_vec;
    else
        obj.exp_info.frame_time = 0;
        obj.exp_info.time_series = 0;
    end
catch
    
    %User defined time_interval. SECONDS!!!!
    obj.exp_info.time_series = [0:obj.exp_info.t_frames-1]*time_interval;
end


end




