 %Go through the raw LSM file and perform the spot tracking in 2D or 3D. 
 %This version uses a LoG filter (same as data-doctor routine) to detect
 %spots. This should work much better than DoG. 

 
 %6-7-18: changing to deal with 3D segmentation.
 
function obj = spot_tracking(obj, params, step)


debug = 1;


Z = obj.exp_info.z_planes;
T = obj.exp_info.t_frames;

%Get segmentation files 
    seg_files = dir([obj.exp_info.nuc_seg_dir,'*.mat']);
    
%LSM time series goes Z first, then time. For each time point, load the
%Z-stack at time t, then max project and segment. 
for t = 1:T

    %Depending on image type, adjust indices
    [~,~,img_type] = fileparts(obj.exp_info.img_file);
    
    switch img_type
        
        case '.lsm'
            %Calculate indices. Skip every other b/c thumbnails? 
            start = Z.*(t-1)*2 + 1;
            indices  = start:2:start + 2*Z-1;
            %Grab this Z-stack
            IMG = tiffread( obj.exp_info.img_file, indices );
            
            %Getting size of example image
            if(iscell(IMG(1).data))
                [y,x] = size( IMG(1).data{1} );
            else
                [y,x] = size( IMG(1).data );
            end

            %Now project
            img = zeros(y,x,Z);
            for z = 1:Z
                if(isfield(params,'seg_channel'))
                    img(:,:,z) = IMG(z).data{params.seg_channel};
                else
                    img(:,:,z) = IMG(z).data;
                end
            end

        case '.tif' %Assuming stitched in image j...
            %Tiff images are z, then channel. T? haven't dealt with these
            %yet, so just ignore t for now. Assume t = 1;
            n_channels = obj.exp_info.Channels;
            start = params.seg_channel;
            indices = start:2:start + 2*Z -1;
            IMG = tiffread( obj.exp_info.img_file, indices);
            %Getting size of example image
            if(iscell(IMG(1).data))
                [y,x] = size( IMG(1).data{1} );
            else
                [y,x] = size( IMG(1).data );
            end

            %Now project
            img = zeros(y,x,Z);
            for z = 1:Z
                img(:,:,z) = IMG(z).data;
            end            
    end
    
    
    %Max project
    mip = max(img,[],3);
    
    
    %LoG filtering (take negative of log filter?)
    img_filter = -im_log_filt_py(mip, params.log_sigma );
    
    %Now make temporary image based on thresholded mip
    im_threshold = imbinarize( img_filter, params.threshold );

    %Collect stats on regions
    stats = regionprops(im_threshold,mip,'centroid','MeanIntensity');
    
    %Check if there were any ROI's
    if( length( stats )==0 )
        continue
    end
        
    
    if(debug)
        figure(7)
        subplot(1,3,1)
        imagesc(im_threshold)
        colormap gray
        %freezeColors
        
        subplot(1,3,2)
        imagesc(img_filter)
        colorbar
        
        
        %Label the regions
        subplot(1,3,1)
        hold on
        for i = 1:length(stats)
            text(stats(i).Centroid(1),stats(i).Centroid(2),num2str(i),'Color','w')
        end
        hold off
        
        subplot(1,3,3)
        imshow( mip,[0,10000]);
        hold on
        viscircles( cat(1,stats.Centroid),4*ones(length(stats),1))
        colormap gray
        pause
    end
    
    
    %Now for this z-stack, load the cell segmentation data. 
    load([obj.exp_info.nuc_seg_dir,seg_files(t).name], 'frame_obj');
    
    %Decide which fitting processes to perform
    switch step.Fit
        
        case '3D'
           
            %Perform fitting in 3D
            fit = fit_this_frame( img, stats, frame_obj, params, step );

            
        case '2D'
    
            %Use maximum projection image or guassial filtered image
            fit = fit_this_frame_mip( mip, stats, frame_obj, params, step );
    end
    
    
    %Append fitting results to frame_obj
    frame_obj.fit = fit;
    save([obj.exp_info.nuc_seg_dir,seg_files(t).name], 'frame_obj','-append');
    
    struct2table(fit) %,'VariableNames',{'Cellid','sigmaXY','sigmaZ','Y','X','Z','Int','BG','D2P'})
    display(['Completed frame: ',num2str(t)])
end



end


