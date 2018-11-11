%Function for making a down-sampled time-lapse image for reviewing cell
%tracking. First max projection, then convert to 8-bit, write tiff
%sequence. 

function obj = create_ds_time_series(obj, params)

debug = 0;


Z = obj.exp_info.z_planes;
T = obj.exp_info.t_frames;


%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = bfGetReader(obj.exp_info.img_file);
series = 1;
BITS = reader.getBitsPerPixel
IMG = obj.exp_info.img_file
%Get the image size of this series. 
size_x = reader.getSizeX;
size_y = reader.getSizeY;

%Output image. 
out_img = zeros(size_y,size_x,T,'uint8' );
[d,fname,ex] = fileparts(obj.exp_info.img_file);

%LSM time series goes Z first, then time. For each time point, load the
%Z-stack at time t, then max project, append to output image. ALways convert to 8-bit. 
disp_str='';
for t = 1:T

    %Get the images for this z-stack according to bioformats reader
    stack = get_stack( reader, t, params.channel);
    
    %Downsample if 16 bit. 
    if(BITS==16)
        out_img(:,:,t) = uint8( max(stack,[],3) ./ 65535 .* 255 );
    elseif(BITS==12)
        out_img(:,:,t) = uint8( max(stack,[],3) ./4095 .*255 ); 
    elseif(BITS==8)
        out_img(:,:,t) = max(stack,[],3);

    else
        disp('unsupported bit type!');
    end 
    
    fprintf(repmat('\b',[1,length(disp_str)+1]))
    disp_str = ['at frame: ',num2str(t)];
    disp(disp_str)
end

%Apparently we don't need to re-order to 5D image with xy? Using order: XYTZC

%temporary save to local disk. 
try
    out_file_name = fullfile('/home/jan/TEMP/',[fname,'_maxp','.tif']);
    %BF exporter. 
    bfsave( out_img, out_file_name,'XYTZC' );
catch
    out_file_name=fullfile('/Users/franklin/Desktop/TEMP/',[fname,'_maxp','.tif']);
    %BF exporter. 
    bfsave( out_img, out_file_name,'XYTZC' );
end

%Final location
final_destination = fullfile(d,[fname,'_maxp','.tif']);

%Now move to final destination. 
movefile(out_file_name, final_destination);
%Add info to exp_info. 
obj.exp_info.max_p_img=final_destination;
obj.save;
end