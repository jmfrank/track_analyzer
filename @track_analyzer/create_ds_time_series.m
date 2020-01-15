%Function for making a down-sampled time-lapse image for reviewing cell
%tracking. First max projection, then convert to 8-bit, write tiff
%sequence. 

function obj = create_ds_time_series(obj, params)

debug = 0;
if nargin<2
    params.channel=1;
    params.z_planes='all';
elseif ~isfield( params, 'z_planes' )
   params.z_planes='all'; 
end

%Figure out how many img_files. 
if iscell(obj.exp_info.img_file)
    
    IMG_files = obj.exp_info.img_file;
    
else
    
    IMG_files{1} = obj.exp_info.img_file;
    
end


%Now loop over files. 
for i = 1:length(IMG_files)
    
    this_file = IMG_files{i};


    %Generate reader. FOR NOW, assume we are looking in series 1. 
    [reader,X,Y,Z,C,T] = bfGetReader(this_file);

    series = 1;
    BITS = reader.getBitsPerPixel;
    IMG = obj.exp_info.img_file;
    %Get the image size of this series. 
    size_x = reader.getSizeX;
    size_y = reader.getSizeY;

    %Output image. 
    out_img = zeros(size_y,size_x,T,'uint8' );
    [d,fname,ex] = fileparts(this_file);

    %LSM time series goes Z first, then time. For each time point, load the
    %Z-stack at time t, then max project, append to output image. ALways convert to 8-bit. 
    disp_str='';
    for t = 1:T

        %Get the images for this z-stack according to bioformats reader
        stack = get_stack( reader, t, params.channel);
        
        % Remove planes as desired. 
        stack = stack(:,:, params.z_planes);
        
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
    %try
    %    out_file_name = fullfile('/home/jan/TEMP/',[fname,'_maxp','.tif']);
    %    %BF exporter. 
    %    bfsave( out_img, out_file_name,'XYTZC' );
    %catch
    %    out_file_name=fullfile('/Users/franklin/Desktop/TEMP/',[fname,'_maxp','.tif']);
    %    %BF exporter. 
    %    bfsave( out_img, out_file_name,'XYTZC' );
    %end

    %Final location
    final_destination = fullfile(d,[fname,'_maxp','.tif']);

    bfsave( out_img, final_destination, 'XYTZC' );
    %Now move to final destination. 
    %movefile(out_file_name, final_destination);
    
    %Out_list. 
    max_p_img_list{i} = final_destination;

end


%Add info to exp_info. 
if i == 1
    obj.exp_info.max_p_img=final_destination;
else
    obj.exp_info.max_p_img = max_p_img_list;
end

obj.save;
end