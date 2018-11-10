 %Go through the raw LSM file and perform the spot tracking in 2D or 3D. 
 %This version uses a LoG filter (same as data-doctor routine) to detect
 %spots. This should work much better than DoG. 

 
 %6-7-18: changing to deal with 3D segmentation.
function varargout = spot_tracking_helper(obj, params, T)

Z = obj.exp_info.z_planes;

%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = bfGetReader(obj.exp_info.img_file);
series = 1;

%Get the image size of this series. 
size_x = reader.getSizeX;
size_y = reader.getSizeY;

%Get segmentation files 
    seg_files = dir([obj.exp_info.nuc_seg_dir,'*.mat']);
    
    
    %Get the images for this z-stack according to bioformats reader
    img = zeros(size_y,size_x,Z);
    for i = 1:Z
        this_plane = reader.getIndex(i-1,params.seg_channel-1,T-1)+1;
        img(:,:,i) = bfGetPlane(reader,this_plane);
    end
    varargout{1} = img;

    if( nargout > 1)
        %New LoG filter in 3D. Using fast implementation.
        img_filter = - log_filter_3D(img, params.log_sigma);   


        %Now open imshow3D for helping to find threshold
        %figure(22)
        %imshow3D_filter(img_filter);
        varargout{2} = img_filter;
    end
    
end


