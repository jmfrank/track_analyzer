%Go through the raw LSM file and perform the spot tracking in 2D or 3D. 
%This revision will detect spots by subtracting the mean intensity value of
%the nuclei using the masks, then look for spots. This ensures less
%influence by variation in nuclei intensities. 

function obj = spot_tracking_r2(obj, params, step)


debug = 1;


Z = obj.exp_info.z_frames;
T = obj.exp_info.t_frames;

%Get filters
    %guassian
    gauss  = fspecial('gaussian',params.pre_gauss);
%Get segmentation files 
    seg_files = dir([obj.exp_info.nuc_seg_dir,'*.tif.mat']);
    
%LSM time series goes Z first, then time. For each time point, load the
%Z-stack at time t, then max project and segment. 
for t = 1:T
    %Calculate indices. Skip every other b/c thumbnails? 
    start = Z.*(t-1)*2 + 1;
    indices  = start:2:start + 2*Z-1;
    %Grab this Z-stack
    IMG = tiffread( obj.exp_info.lsm_file, indices );
    
    [y,x] = size( IMG(1).data );
    %Now collect z-stack 
    img = zeros(y,x,Z);
    for z = 1:Z
        img(:,:,z) = IMG(z).data;
    end

    %Project this z-stack. 
    mip = max(img,[],3);
    %Smooth / Filter / get blobs
    gauss_diff = mip - imgaussfilt(mip, 5);
    gauss_final = imgaussfilt(gauss_diff,1);
    
    im_threshold = imbinarize( gauss_final, params.threshold );
    


    %Now make temporary image based on thresholded mip
    stats = regionprops(im_threshold,mip,'centroid','MeanIntensity');
    
    if(debug)
        figure(7)
        subplot(1,2,1)
        imagesc(im_threshold)
        colormap gray

        subplot(1,2,2)
        imagesc(gauss_final)
        colormap gray
        %Label the regions
        subplot(1,2,1)
        hold on
        for i = 1:length(stats)
            text(stats(i).Centroid(1),stats(i).Centroid(2),num2str(i),'Color','w')
        end
        hold off
        pause
    end
    

    %Now for this z-stack, load the cell segmentation data. 
    load([obj.exp_info.nuc_seg_dir,seg_files(t).name], 'frame_obj');
    
    %binary image of this frame
    bin = frame_obj.BW;
    %Get the mean intensity of each cell
    stats_bg = regionprops(bin,mip,'MeanIntensity','PixelIdxList');
    %Now create background image
    bg = zeros(y,x);
    for i = 1:length(stats_bg)
        this_list = stats_bg(i).PixelIdxList;
        this_int  = stats_bg(i).MeanIntensity;
        bg(this_list) = this_int;
    end
    %Now subtract background from mip
    mip_corrected = mip - bg;
    
    if(debug)
    figure(8)
    subplot(1,2,1)
    imagesc(bin)
    colormap gray
    subplot(1,2,2)
    imagesc(mip_corrected);
    colormap gray
    pause
    end
        
    
           
    %Perform fitting in 3D
    %results = fit_this_frame( img, stats, frame_obj, params, step );
    
    %Use maximum projection image or guassial filtered image
    results = fit_this_frame_mip( gauss_final, stats, frame_obj, params, step );
    
    %Append results to frame_obj
    frame_obj.results = results;
    save([obj.exp_info.nuc_seg_dir,seg_files(t).name], 'frame_obj','-append');
    
    array2table(results,'VariableNames',{'Cellid','sigmaXY','sigmaZ','Y','X','Z','Int','BG','D2P'})
end



end

function results = fit_this_frame_mip( mip, stats, frame_obj, params, step )

%Get dimensions of stack
[l,w] = size(mip); 

%Now loop through cells and determine if there's a spot in each cell at
%this time
n_cells = length(frame_obj.refined_cells);
spot_count = 0;

centers = cat(1,stats.Centroid);
int_vals = cat(1,stats.MeanIntensity);

for i = 1:n_cells

    %Load contour 
    cn = frame_obj.refined_cells{i};
    %Expand contour
       x = cn(:,1);
       y = cn(:,2);

       %Transform
       A = params.expand_ratio/100 + 1;
       x = A*x+(1-A)*mean(x);
       y = A*y+(1-A)*mean(y);
       
       %Replace
       cn = [x,y];
    %Create binary
    
    %Search which blobs are within this cell MASK
    IN = inpolygon(centers(:,1),centers(:,2),cn(:,1),cn(:,2));
    %Count members
    n_spots = sum(IN);

    %If there are any spots in this cell, analyze / save
    if( n_spots > 0 ) 
        %Let's choose the maximum valued spot
        
            %Get the index
            idx = find(IN==1);
            %Intensity values
            iv = int_vals(idx);
            
            %Get highest spot id
            [~,jdx] = max( iv );
            idx = idx(jdx);
            
            spot_count = spot_count + 1;
            ctr = round(centers(idx,:));

            %Create bounding box
            x = ctr(1)-params.search_radius:ctr(1)+params.search_radius;
            sel = x > 0 & x <= l;
            x = x(sel);
            y = ctr(2)-params.search_radius:ctr(2)+params.search_radius;
            sel = y > 0 & y <= w;
            y = y(sel);

            %Select the roi
            search_stack = mip(y,x);

            %Guess center for fitting
            ctr_guess = [params.search_radius+1,params.search_radius+1];        

            %Fixed sigma
            if(step.FixedSigma)
                %Fit search stack to 2D gaussian
                vals = fit_2D_gauss_fixed_sigma(search_stack,ctr_guess,params.gauss);
                %Manually add back sigma values
                vals = [ params.gauss.xy_sigma_est,vals];
            else
                %Fit search stack to 2D gaussian
                vals = fit_2D_gauss(search_stack,ctr_guess,params.gauss);
            end
            
            %
            fit(spot_count,:) = [vals(1),0,vals(2:3),0,vals(4:5)];
            
            %Calculate distance to nuclear periphery
            spot_loc = fit(spot_count,4:-1:3);
            dist_vec = dist2(spot_loc,cn).^0.5;
            [min_dist(spot_count), min_idx] = min(dist_vec);

            %Keep track of which cell has this spot
            cell_id(spot_count) = i;
            
%             figure(10);
%             title(['Spot ID: ',num2str(idx)])
%             hold on
%             mip_tmp = max(search_stack,[],3);
%             imagesc(mip_tmp)
%             plot(fit(spot_count,4),fit(spot_count,3),'*r')
%             hold off
%             set(gca,'ydir','reverse')
%             colormap gray
%             pause
% %   
%             figure(8);
%             imshow3D(search_stack)
%             hold on
%             plot(fit(spot_count,4),fit(spot_count,3),'*r')
%             hold off
% 
%             fit(spot_count,3) = fit(spot_count,3) + y(1) -1;
%             fit(spot_count,4) = fit(spot_count,4) + x(1) -1;            
%             array2table(fit(spot_count,:),'VariableNames',{'xy_sigma','z_sigma','Y','X','Z','Int','BG'})
%             
%             %Shift x-y vals back into frame reference space
% 
%             pause
    end       
end

%Compile data [ cell id, fitting measurements, dist to periphery]
results = [cell_id',fit, min_dist'];

end

function results = fit_this_frame( img, stats, frame_obj, params, step )

%Get dimensions of stack
[l,w,h] = size(img); 

%Now loop through cells and determine if there's a spot in each cell at
%this time
n_cells = length(frame_obj.refined_cells);
spot_count = 0;

centers = cat(1,stats.Centroid);
int_vals = cat(1,stats.MeanIntensity);

for i = 1:n_cells

    %Load contour 
    cn = frame_obj.refined_cells{i};
    %Expand contour
       x = cn(:,1);
       y = cn(:,2);

       %Transform
       A = params.expand_ratio/100 + 1;
       x = A*x+(1-A)*mean(x);
       y = A*y+(1-A)*mean(y);
       
       %Replace
       cn = [x,y];
    %Create binary
    
    %Search which blobs are within this cell MASK
    IN = inpolygon(centers(:,1),centers(:,2),cn(:,1),cn(:,2));
    %Count members
    n_spots = sum(IN);

    %If there are any spots in this cell, analyze / save
    if( n_spots > 0 ) 
        %Let's choose the maximum valued spot
        
            %Get the index
            idx = find(IN==1);
            %Intensity values
            iv = int_vals(idx);
            
            %Get highest spot id
            [~,jdx] = max( iv );
            idx = idx(jdx);
            
            spot_count = spot_count + 1;
            ctr = round(centers(idx,:));

            %Create bounding box
            x = ctr(1)-params.search_radius:ctr(1)+params.search_radius;
            sel = x > 0 & x <= l;
            x = x(sel);
            y = ctr(2)-params.search_radius:ctr(2)+params.search_radius;
            sel = y > 0 & y <= w;
            y = y(sel);

            %Select the roi (use all z)
            search_stack = img(y,x,:);

            %Guess z-center using location of the maximum pixel
            [max_val,max_idx] = max(search_stack(:));
            [~,~,ctr_z] = ind2sub(size(search_stack),max_idx);
            %[~,ctr_z] = max(sm_img(ctr(2),ctr(1),:));
            z = ctr_z; %ctr_z - params.search_z: ctr_z + params.search_z;
            sel = z > 0 & z <= h;
            z = z(sel);

            %Guess center for fitting
            ctr_guess = [params.search_radius+1,params.search_radius+1,z];        

            %Fixed sigma
            if(step.FixedSigma)
                %Fit search stack to 3D gaussian
                vals = fit_3D_gauss_fixed_sigma(search_stack,ctr_guess,params.gauss);
                %Manually add back sigma values
                fit(spot_count,:) = [ params.gauss.xy_sigma_est,params.gauss.z_sigma_est,vals];
            else
                %Fit search stack to 3D gaussian
                fit(spot_count,:) = fit_3D_gauss(search_stack,ctr_guess,params.gauss);
            end

            %Calculate distance to nuclear periphery
            spot_loc = fit(spot_count,4:-1:3);
            dist_vec = dist2(spot_loc,cn).^0.5;
            [min_dist(spot_count), min_idx] = min(dist_vec);

            %Keep track of which cell has this spot
            cell_id(spot_count) = i;
            
%             figure(10);
%             title(['Spot ID: ',num2str(idx)])
%             hold on
%             mip_tmp = max(search_stack,[],3);
%             imagesc(mip_tmp)
%             plot(fit(spot_count,4),fit(spot_count,3),'*r')
%             hold off
%             set(gca,'ydir','reverse')
%             colormap gray
% 
%             figure(8);
%             imshow3D(search_stack)
%             hold on
%             plot(fit(spot_count,4),fit(spot_count,3),'*r')
%             hold off
% 
%             fit(spot_count,3) = fit(spot_count,3) + y(1) -1;
%             fit(spot_count,4) = fit(spot_count,4) + x(1) -1;            
%             array2table(fit(spot_count,:),'VariableNames',{'xy_sigma','z_sigma','Y','X','Z','Int','BG'})
%             
%             %Shift x-y vals back into frame reference space
% 
%             pause
    end       
end

%Compile data [ cell id, fitting measurements, dist to periphery]
results = [cell_id',fit, min_dist'];

end

