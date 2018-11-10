 %Go through the raw LSM file and perform the spot tracking in 2D or 3D. 
 %This version uses a LoG filter (same as data-doctor routine) to detect
 %spots. This should work much better than DoG. 

 
 %6-7-18: changing to deal with 3D segmentation.
 
function obj = spot_tracking(obj, params, step)

Z = obj.exp_info.z_planes;
T = obj.exp_info.t_frames;


%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = bfGetReader(obj.exp_info.img_file);
series = 1;

%Get the image size of this series. 
size_x = reader.getSizeX;
size_y = reader.getSizeY;

%Get segmentation files 
seg_files = obj.get_frame_files;

%LSM time series goes Z first, then time. For each time point, load the
%Z-stack at time t, then max project and segment. 
for t = 1:T

    
    %Get the images for this z-stack according to bioformats reader
    img = zeros(size_y,size_x,Z);
    for i = 1:Z
        this_plane = reader.getIndex(i-1,params.seg_channel-1,t-1)+1;
        img(:,:,i) = bfGetPlane(reader,this_plane);
    end
    
    %IMPORTANT: we often deal with stitched images. So we should fill in
    %zero-valued pixels. 
    img = fill_edges(img);
    
    %New LoG filter in 3D. Using fast implementation.
    img_filter = - log_filter_3D(img, params.log_sigma);

    %Binarize image based on threshold. 
    im_threshold = img_filter >= params.threshold;

    %Load cell data and collect centroids. 
    load(seg_files{t}, 'frame_obj');
    
    %Now just multily by the cell nucleus mask
        %Check if BW is 2D or 3D. 
        if(numel(size(frame_obj.BW))==2)
            BW = repmat(frame_obj.BW,[1,1,size(im_threshold,3)]);
        else
            BW = frame_obj.BW;
        end
    im_threshold = im_threshold.*BW;
        
    %Collect stats on regions
    stats = regionprops(logical(im_threshold),'Centroid');
    

    %Estimate assignment. 
    cell_centroids = cat(1,frame_obj.centroids{:});
    spot_centroids = cat(1,stats.Centroid);
    if(size(cell_centroids,2)==2)
        spot_centroids = spot_centroids(:,1:2);
    end
           
    D = pdist2(spot_centroids, cell_centroids);
    [~,assignment] = min(D,[],2);
    A = num2cell(assignment);
    [stats.assignment] = A{:};
    
    %Refine assigment for those that are too close. Compare distance to nearest and next-nearest. If within a threshold,
        %refine the assignment. 
        [D_sort,sort_I] = sort(D,2);
        cmp_dist = D_sort(:,2)- D_sort(:,1);
        these_conflict = find(cmp_dist < 15); %15 pixels min distance to next neighbor. 
        cell_ids       = sort_I(these_conflict,1:2);
        %Empty reference image. 
        E = zeros(size(BW,1),size(BW,2));

        for c = 1:length(these_conflict)
           %this spot centroid. 
           this_spot_centroid = spot_centroids(these_conflict(c),:);

           %Calculate contours for the two possible cells. 
           these_cells = cell_ids(c,:);

           for z = 1:length(these_cells)
                cell_img = E;    
                px = frame_obj.PixelIdxList{ these_cells(z)};
                cell_img(px) = 1;
                [x,y] = C2xyz(contourc(cell_img,[0.5,0.5]));
                IN(z) = inpolygon( this_spot_centroid(1),this_spot_centroid(2),x{1},y{1});
           end

           %Check if spot is inside at least one cell. 
           if(sum(IN)==1)
               %Assign spot to correct cell. 
               stats(these_conflict(c)).assignment = these_cells(IN);
           else
               disp('Did not assign spot')
           end

        end


    %Check if there were any empty ROI's
    if( length( stats )==0 )
        continue
    end
    
    if(step.debug)
        figure(7)
        imshow3D(img_filter)
        colorbar
        axis equal
        axis tight
        hold on
        %Concat xy coordinates of centroids
        centroids = cat(1,stats.Centroid);
        centroids = centroids(:,[1,2]);
        viscircles( centroids,4*ones(length(stats),1))
        colormap gray
        pause
    end
    
    
    
    %Decide which fitting processes to perform
    switch step.Fit
        
        %3d gaussian. 
        case '3D'
           
            %Perform fitting in 3D
            fit = fit_this_frame( img, stats, params, step );
           
        %2d gaussian
        case '2D'
            
            %Max project
            mip = max(img,[],3);
            %Use maximum projection image or guassial filtered image
            fit = fit_this_frame_mip( mip, stats, frame_obj, params, step );
    
        case 'Centroid'
            
            %Simple centroid fitting scheme
            fit = fit_this_frame_centroid( img, stats, params);

    end
    
    %Append fitting results to frame_obj
    frame_obj.fit = fit;
    save(seg_files{t}, 'frame_obj','-append');
    
    struct2table(fit) %,'VariableNames',{'Cellid','sigmaXY','sigmaZ','Y','X','Z','Int','BG','D2P'})
    display(['Completed frame: ',num2str(t)])
end



end


