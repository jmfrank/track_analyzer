 %Spot tracking without caring about cells.  
function obj = spot_tracking_no_cells(obj, params, step)

Z = obj.exp_info.z_planes;
T = obj.exp_info.t_frames;

%Loop over multiple files (pre-bleach and post-bleach). 
files = {obj.exp_info.pre_bleach, obj.exp_info.img_file};
seg_files_all = {};
for F = 1:length(files)
    
   [obj,seg_files] = inner_function(obj, F, files{F}, params, step);
    seg_files_all = [seg_files_all, seg_files];
end

%Now add list of seg_files to exp_info...
obj.exp_info.seg_files = seg_files_all;
end

%% inner function for running pipeline. 
function [obj, seg_files] = inner_function( obj, f_id, img_file,params,step)
        
%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = bfGetReader(img_file);
series = 1;

%Get the image size of this series. 
size_x = reader.getSizeX;
size_y = reader.getSizeY;
T = reader.getSizeT;
Z = reader.getSizeZ;

%LSM time series goes Z first, then time. For each time point, load the
%Z-stack at time t, then max project and segment. 
seg_files = {};
for t = 1:T
    
    %Get the images for this z-stack according to bioformats reader
    img = zeros(size_y,size_x,Z);
    for i = 1:Z
        this_plane = reader.getIndex(i-1,params.seg_channel-1,t-1)+1;
        img(:,:,i) = bfGetPlane(reader,this_plane);
    end
    
    %Change process for 2D or 3D
    switch step.Fit
        
        case '3D'
            %New LoG filter in 3D. Using fast implementation.
            img_filter = -log_filter_3D(img, params.log_sigma);

            %Try smoothing image
            %img_filter = imgaussfilt3(img, params.gauss_sigma);

            img_bw = img_filter >= params.thresh;

            stats = regionprops(logical(img_bw),'Centroid','Area');
            %Filter out ridiculous values. 
            sel = [stats.Area] < 30;
            stats = stats(sel);
            A = num2cell( ones(length(stats),1));
            [stats.assignment] = A{:};
            
            %Perform fitting in 3D
            fit = fit_this_frame( img, stats, params, step );

        case '2D'
            
            mip = max(img,[],3);
            h = fspecial('log',[15,15],8); %,[5,5],params.log_sigma(1));
            img_filter = -imfilter(mip,h,'symmetric', 'conv');
            
            img_bw = img_filter >= params.thresh;
            stats = regionprops(img_bw,'Centroid','Area');
            A = num2cell( ones(length(stats),1));
            [stats.assignment] = A{:};
            %Use maximum projection image or guassial filtered image
            fit = fit_this_frame_mip( mip, stats, params, step );
            
        case 'Centroid'
            
            mip = max(img,[],3);
            h = fspecial('log',[15,15],8); %,[5,5],params.log_sigma(1));
            img_filter = -imfilter(mip,h,'symmetric', 'conv');
            
            img_bw = img_filter >= params.thresh;
            stats = regionprops(img_bw,'Centroid','Area');
            A = num2cell( ones(length(stats),1));
            [stats.assignment] = A{:};    
            if(isempty(stats))
                fit = [];
            else
                fit = fit_this_frame_centroid_2D( mip, stats, params );
            end
    end
       
    if(step.debug)
        figure(7)
        clf
        imshow3D_spot_filter(img_filter,[-0.1,0.3],params.thresh);
        hold on
        %Also plot the ROI
        %p = [obj.exp_info.ROI.x_pos,obj.exp_info.ROI.y_pos];
        %viscircles(p,50,'color','g')
        %if(length(stats)>0)
        %    %Concat xy coordinates of centroids
        %    centroids = cat(1,stats.Centroid);
        %    centroids = centroids(:,[1,2]);
        %    viscircles( centroids,8*ones(length(stats),1))
        %    colormap gray
        %end
        pause
    end
    
    %Saving 
    seg_files{t} = [obj.exp_info.nuc_seg_dir,'file_',num2str(f_id),'_frame_',sprintf('%04d',t),'.mat'];
    
    %Check if there were any ROI's
    if( length( stats )==0 )
        %Save empty fit structure
        frame_obj.fit = [];
    else
        %Append fitting results to frame_obj
        frame_obj.fit = fit;
    end
    
    if(exist(seg_files{t},'file'))
        save(seg_files{t}, 'frame_obj','-append');
    else
        save(seg_files{t},'frame_obj')
    end
    
    %struct2table(fit) %,'VariableNames',{'Cellid','sigmaXY','sigmaZ','Y','X','Z','Int','BG','D2P'})
    display(['Completed frame: ',num2str(t)])
end

reader.close();
end


