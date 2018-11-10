 %This function does not include cell segmentation step. Useful for single
 %cell images (i.e. FRAP of nascent mRNA spots)

function obj = spot_tracking_only(obj, params, step)


debug = 1;


Z = obj.exp_info.z_frames;
T = obj.exp_info.t_frames;

%We will write directory to the obj.results cause we don't need to worry
%about keep track of cell segmentation. Cell_id = 1 for all localizations. 
obj.results = [];map

%Get filters
    %guassian
    gauss  = fspecial('gaussian',params.pre_gauss);
%Get segmentation files 
    seg_files = dir([obj.exp_info.nuc_seg_dir,'*.mat']);
    
%LSM time series goes Z first, then time. For each time point, load the
%Z-stack at time t, then max project and segment. 
for t = 31 %1:T
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
        imagesc(mip)
        colormap gray
        pause
    end
    

    %Figure out if we're using 3D or 2D
    switch step.fit_type
        
        case '3D'            
            %Perform fitting in 3D
            results = fit_this_frame( img, stats, frame_obj, params, step );
        case '2D'
            %Use maximum projection image or guassial filtered image
            fit = fit_this_frame_mip( mip, stats, frame_obj, params, step );
    end
    
    %Append results to frame_obj
    frame_obj.fit = fit;
    save([obj.exp_info.nuc_seg_dir,seg_files(t).name], 'frame_obj','-append');
    
    struct2table(fit) %,'VariableNames',{'Cellid','sigmaXY','sigmaZ','Y','X','Z','Int','BG','D2P'})
    display(['Completed frame: ',num2str(t)])
end



end



%Pass a 2D frame, and find spots
function fit = fit_this_frame_mip( mip, stats, frame_obj, params, step )

debug = 0;
%Get dimensions of stack
[l,w] = size(mip); 

%Now loop through cells and determine if there's a spot in each cell at
%this time
n_cells = length(frame_obj.refined_cells);
spot_count = 0;

%Get vectors of peak centers and mean-intensities
centers = cat(1,stats.Centroid);
int_vals = cat(1,stats.MeanIntensity);

%Empty structure for fitting measurements. Input is the number of
%centroids. Will need to remove empty entries later on. 
fit = gen_fit_struct( length( centers) ); 

%Loop over cells
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
                %Fit search stack to 2D gaussian (this fit doesn't exist
                %yet?!
                fit(spot_count) = fit_2D_gauss_fixed_sigma(search_stack,ctr_guess,params.gauss);
 
            else
                %Fit search stack to 2D gaussian (independent sigma x,y)
                fit(spot_count) = fit_2D_gauss_2sigmas(search_stack,ctr_guess,params.gauss);
            end
            
            
            %Shift x-y vals back into frame reference space
            fit(spot_count).y_pos = fit(spot_count).y_pos + y(1) -1;
            fit(spot_count).x_pos = fit(spot_count).x_pos + x(1) -1;
            
            %Calculate distance to nuclear periphery
            spot_loc = [fit(spot_count).x_pos, fit(spot_count).y_pos];
            dist_vec = dist2(spot_loc,cn).^0.5;
            fit(spot_count).dist_2_periphery = min(dist_vec);
            
            %Assign this fit to cell 'i'
            fit(spot_count).cell_id = i;
            
            if(debug)
                figure(10);
                title(['Spot ID: ',num2str(idx)])
                hold on
                mip_tmp = max(search_stack,[],3);
                imagesc(mip_tmp)
                plot(fit(spot_count,4),fit(spot_count,3),'*r')
                hold off
                set(gca,'ydir','reverse')
                colormap gray
                pause
                figure(8);
                imshow3D(search_stack)
                hold on
                plot(fit(spot_count,4),fit(spot_count,3),'*r')
                hold off
                array2table(fit(spot_count,:),'VariableNames',{'xy_sigma','z_sigma','Y','X','Z','Int','BG'})
                pause
            end
    end       
end

%Now we remove empty entries of the fit structure
fit = fit(1:spot_count);

end

%Pass a 3D frame and find spots
function results = fit_this_frame( img, stats, frame_obj, params, step )

%Get dimensions of stack
[l,w,h] = size(img); 

%Now loop through cells and determine if there's a spot in each cell at
%this time
n_cells = length(frame_obj.refined_cells);
spot_count = 0;

%Get vectors of peak centers and mean-intensities
centers = cat(1,stats.Centroid);
int_vals = cat(1,stats.MeanIntensity);

%Empty structure for fitting measurements. Input is the number of
%centroids. Will need to remove empty entries later on. 
fit = gen_fit_struct( length( centers) ); 

%Loop over cells
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
            z = ctr_z; 
            sel = z > 0 & z <= h;
            z = z(sel);

            %Guess center for fitting
            ctr_guess = [params.search_radius+1,params.search_radius+1,z];        

            %Fixed sigma
            if(step.FixedSigma)
                %Fit search stack to 3D gaussian
                fit(spot_count) = fit_3D_gauss_fixed_sigma(search_stack,ctr_guess,params.gauss);
            else
                %Fit search stack to 3D gaussian
                fit(spot_count) = fit_3D_gauss(search_stack,ctr_guess,params.gauss);
            end
            
            %Shift x-y vals back into frame reference space
            fit(spot_count).y_pos = fit(spot_count).y_pos + y(1) -1;
            fit(spot_count).x_pos = fit(spot_count).x_pos + x(1) -1;
            
            %Calculate distance to nuclear periphery
            spot_loc = [fit(spot_count).x_pos, fit(spot_count).y_pos];
            dist_vec = dist2(spot_loc,cn).^0.5;
            fit(spot_count).dist_2_periphery = min(dist_vec);
            
            %Assign this fit to cell 'i'
            fit(spot_count).cell_id = i;
            
            if(debug)
                figure(10);
                title(['Spot ID: ',num2str(idx)])
                hold on
                mip_tmp = max(search_stack,[],3);
                imagesc(mip_tmp)
                plot(fit(spot_count,4),fit(spot_count,3),'*r')
                hold off
                set(gca,'ydir','reverse')
                colormap gray

                figure(8);
                imshow3D(search_stack)
                hold on
                plot(fit(spot_count,4),fit(spot_count,3),'*r')
                hold off

                array2table(fit(spot_count,:),'VariableNames',{'xy_sigma','z_sigma','Y','X','Z','Int','BG'})

                pause
            end
    end       
end

%Now we remove empty entries of the fit structure
fit = fit(1:spot_count);

end


