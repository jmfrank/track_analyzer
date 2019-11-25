 %Go through the raw LSM file and perform the spot tracking in 2D or 3D. 
 %This version uses a LoG filter (same as data-doctor routine) to detect
 %spots. This should work much better than DoG. 

 
 %6-7-18: changing to deal with 3D segmentation.
 %6-29-18: developed a way to locally threshold cells based using
 %percentile of intensities for each nucleus. 
 
function obj = spot_tracking_LIVE(obj, params, step, frames)

% Image channel
channel_str = ['seg_channel_',pad(num2str(params.seg_channel),2,'left','0')];

%Generate reader. FOR NOW, assume we are looking in series 1. 
[reader,X,Y,Z,C,T] = bfGetReader(obj.exp_info.img_file);
series = 1;

%Figure out if call to function included specific frames. 
if nargin < 4
    frames = 1:T;
else
    frames = frames;
end

%Get the image size of this series. 
size_x = reader.getSizeX;
size_y = reader.getSizeY;

%Get segmentation files 
seg_files = dir([obj.exp_info.nuc_seg_dir,'*.mat']);
    
%LSM time series goes Z first, then time. For each time point, load the
%Z-stack at time t, then max project and segment. 
for t = frames
    %Get the images for this z-stack according to bioformats reader
    img = zeros(size_y,size_x,Z);
    for i = 1:Z
        this_plane = reader.getIndex(i-1,params.seg_channel-1,t-1)+1;
        img(:,:,i) = bfGetPlane(reader,this_plane);
    end
    
    %Load frame obj to get segmented cells. 
    load([obj.exp_info.nuc_seg_dir,seg_files(t).name], 'frame_obj');
    
    %Optional smoothing step.
    if step.smooth_img
       
        img = imgaussfilt3(img, params.smoothing_sigma);
    end
    
    %% Subtract local background? 
    if step.subtract_local_background 
        
        % Using the nuclear masks, subtract the median/mean intensity value
        % from each cell nuclei. This will prevent large gradients at
        % nuclear edge. 
        img = subtract_local_backgroud(img, frame_obj.(channel_str)) ;
    end
    
    %New LoG filter in 3D. Using fast implementation.
    img_filter = -log_filter_3D(img, params.log_sigma);

    %Load frame obj to get segmented cells. 
    load([obj.exp_info.nuc_seg_dir,seg_files(t).name], 'frame_obj');
    
    %Adding a local thresholding algorithm. Stats contains centroids and cell assignments.  
    stats = local_threshold_cells( img, img_filter, frame_obj, params);
        
    %Check if there were any ROI's
    if( length( stats )==0 )
        continue
    end
    
    % Filter out smaller regions. fit_this_frame_centroid
    if step.max_size
        areas = [stats.Area];
        sel = areas < step.max_size;
        stats = stats(sel);
    end
    
    if step.debug
        figure(10)
        clf
        imshow3D(img_filter);
        hold on
        %Concat xy coordinates of centroids
        centroids = cat(1,stats.Centroid);
        centroids = centroids(:,[1,2]);
        for i = 1:length(stats)
            h = viscircles( centroids(i,:), 15);
            h.Children(1).UserData = i;           
        end
        
        cell_centroids = cat(1,frame_obj.(channel_str).centroids);
        for i = 1:length(cell_centroids)
            text(cell_centroids{i}(1),cell_centroids{i}(2),num2str(i),'fontsize',20,'Color','w')
        end
        
        colormap gray
    end
    
    %Decide which fitting processes to perform
    switch step.Fit
        
        case '3D'
           
            %Perform fitting in 3D
            fit = fit_this_frame( img, stats, params, step );

            
        case '2D'
    
            %Use maximum projection image or guassial filtered image
            fit = fit_this_frame_mip( mip, stats, frame_obj, params, step );
            
        case 'Centroid'
            
            fit = fit_this_frame_centroid( img, stats, params);
            
    end
    
    %% Remove fits if necessary.
    if step.max_fits_per_nucleus
        new_fits = [];
        cell_ids = [fit.cell_id];
        
        for i = unique(cell_ids)
            these_fits = cell_ids == i;
            ints = [fits(these_fits).sum_int];

            % Find by top intensity?
            I = [fit.sum_int];
            [Is,idx] = sort(I,'descend');
            
            if length(Is) > params.max_fits
                new_fits=[new_fits,fits(idx(1:params.max_fits))];
            end
            
        end
    end
    
    
    %Append fitting results to frame_obj
    frame_obj.(channel_str).fit = fit;
    save([obj.exp_info.nuc_seg_dir,seg_files(t).name], 'frame_obj','-append');
    
    %struct2table(fit) %,'VariableNames',{'Cellid','sigmaXY','sigmaZ','Y','X','Z','Int','BG','D2P'})
    display(['Completed frame: ',num2str(t)])
end

reader.close();


end


%% 

function img = subtract_local_backgroud(img, frame_obj)



    % Loop over cell nuclei in BW. 
    n = length(frame_obj.PixelIdxList);
    
    for i = 1:n
        
        bw = false(size(frame_obj.BW));
        
        bw( frame_obj.PixelIdxList{i} ) = 1;
        
        mean_val = mean( img(bw));
        %median_val = median( img(bw));
        
        img(bw) = img(bw)-mean_val;
        
    end
    
    % Set the background to zero?
    bg = ~frame_obj.BW;
    img(bg) = 0;
    
end
        

    

