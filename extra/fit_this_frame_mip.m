%Pass a 2D frame, and find spots. Used for nascent RNA spot localization.

%2-23-18: Stitched images often contain an edge where there's zero
%intensity outside of the actual image. This creates false peak detections
%due to the sharp edge. I added some lines to mask out the peaks too close
%to the such edges. 
%7-12-18: getting rid of cell-based stuff

function fit = fit_this_frame_mip( mip, stats, params, step )

debug = 0;
%Get dimensions of stack
[l,w] = size(mip);

%Get vectors of peak centers and mean-intensities
centers = cat(1,stats.Centroid);

%Make a special shell-mask for querying the local background. 
    bounds=params.gauss.bg_mask;
    B = strel('disk',bounds(2),0);
    A = strel('disk',bounds(1),0);
    
    fgMASK = padarray(A.Neighborhood,(params.search_radius-bounds(1))*[1,1]);
    bgMASK = padarray(B.Neighborhood,(params.search_radius-bounds(2))*[1,1])...
        & ~fgMASK;

%Empty fit structure    
fit = gen_fit_struct(size(centers,1));

for i = 1:size(centers,1)
      disp(['Fitting spot: ',num2str(i)])

    ctr = round(centers(i,:));

    %Create bounding box
    x = ctr(1)-params.search_radius:ctr(1)+params.search_radius;
    sel_x = x > 0 & x <= w;
    x = x(sel_x);
    y = ctr(2)-params.search_radius:ctr(2)+params.search_radius;
    sel_y = y > 0 & y <= l;
    y = y(sel_y);
    %Create fg and bg masks. 
    this_fg = fgMASK(sel_y,sel_x);
    this_bg = bgMASK(sel_y,sel_x);
    
    %Select the roi
    search_stack = double(mip(y,x));
    
    %Intensity guess
    peak = mean(search_stack( this_fg ));
    
    %Background guess. Should do this for each fit cause edge effects :(
    bg_guess = mean( search_stack( this_bg ));
    
    %Provide a background and amplitude guess 
    params.gauss.bg_guess  = double(bg_guess);
    params.gauss.int_guess = double(peak - bg_guess);
    
    %Guess center for fitting
    ctr_guess = [params.search_radius+1,params.search_radius+1];        

    %Fixed sigma
    switch step.fit_type
        case '2D_gauss'
            %Fit search stack to 2D gaussian (independent sigma x,y)
            fit(i) = fit_2D_gauss(search_stack,ctr_guess,params.gauss,this_fg);
            
            %fit(i) = fit_2D_gauss_2sigmas(search_stack,ctr_guess,params.gauss);
        case '2D_blob'
            
            %Just need the intensity above a background value.
            %Estimate background by sampling away from ROI
            fit(i) = fit_2D_blob(search_stack,ctr_guess,params.gauss);
    end

    %Shift x-y vals back into frame reference space
    fit(i).pos = fit(i).pos + [y(1)-1,x(1)-1];


    %Calculate distance to nuclear periphery
    %spot_loc = fit(i).pos([2,1]);
    %dist_vec = dist2(spot_loc,cn).^0.5;
    %fit(i).dist_2_periphery = min(dist_vec);

    %Assign this fit to cell 'i'
    %fit(i).cell_id = 1;

    if(debug)
        figure(10);
        title(['Spot ID: ',num2str(1)])
        hold on
        imagesc(mip)
        plot(fit(i).pos(2),fit(i).pos(1),'*r')
        hold off
        set(gca,'ydir','reverse')
        colormap gray
        pause
    end
    
end
    
end