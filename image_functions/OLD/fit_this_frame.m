%Pass a 3D frame and find spots. Works with 2D or 3D segmentation

% 6-20-18: Improving speed in assigning each spot to a nucleus. Only
% searching if spot is already close to a nucleus. The threshold for
% 'closeness' will be determined by the max volume or area of all cells.
% Assume an ellipse/ellipseoid with 3:1 length ratio. Use max dimension of
% this ellipse to define distance threshold. 
% 6-29-18: changing based on the local threshold step just before this. No
% need to verify if spot is within cell. 
function fig = fit_this_frame( img, stats, params, step )

debug = 0;
%Get dimensions of stack
[l,w,h] = size(img); 

%Get vectors of peak centers and mean-intensities
centers = cat(1,stats.Centroid);

%Empty structure for fitting measurements. Input is the number of
%centroids. Will need to remove empty entries later on. 
fit = gen_fit_struct( size(centers,1) );

%Make a special shell-mask for querying the local background. 
    bounds=params.gauss.bg_mask;
    B = strel('sphere',bounds(2));
    A = strel('sphere',bounds(1));
    
    fgMASK = padarray(A.Neighborhood,(params.search_radius-bounds(1))*[1,1,1]);
    bgMASK = padarray(B.Neighborhood,(params.search_radius-bounds(2))*[1,1,1])...
        & ~fgMASK;

for i = 1:size(centers,1)
    disp(['Fitting spot: ',num2str(i)])
    
    %approximate center. 
    ctr = round(centers(i,:));

    %Create bounding box
    x = ctr(1)-params.search_radius:ctr(1)+params.search_radius;
    sel_x = x > 0 & x <= w;
    x = x(sel_x);
    y = ctr(2)-params.search_radius:ctr(2)+params.search_radius;
    sel_y = y > 0 & y <= l;
    y = y(sel_y);
    z = ctr(3)-params.search_radius:ctr(3)+params.search_radius;
    sel_z = z > 0 & z <= h;
    z = z(sel_z);
    %Select the roi (use all z)
    search_stack = double(img(y,x,z));

    %Guess center for fitting
    ctr_guess = [params.search_radius+1,params.search_radius+1,params.search_radius+1];        
    
    this_fg = fgMASK(sel_y,sel_x,sel_z);
    this_bg = bgMASK(sel_y,sel_x,sel_z);
    
    %Intensity guess
    peak = mean(search_stack( this_fg ));
    
    %Background guess. Should do this for each fit cause edge effects :(
    bg_guess = mean( search_stack( this_bg ));
    
    %Provide a background and amplitude guess 
    params.gauss.bg_guess  = double(bg_guess);
    params.gauss.int_guess = double(peak - bg_guess);
    
    %Fixed sigma
    if(step.FixedSigma)
        %Fit search stack to 3D gaussian
        fit(i) = fit_3D_gauss_fixed_sigma(search_stack,ctr_guess,params.gauss,this_fg);
    else
        %Fit search stack to 3D gaussian
        fit(i) = fit_3D_gauss_SIGxy_constant( search_stack, ctr_guess,params.gauss,this_fg);
    end

    %Shift x-y vals back into frame reference space
    fit(i).pos = fit(i).pos + [y(1)-1,x(1)-1,z(1)-1];
        
    %Assign fit to cell
    fit(i).cell_id = stats(i).assignment;
    
    
    if(step.NuclearSegmentation=='2D')
        %Calculate distance to nuclear periphery
        spot_loc = [fit(i).x_pos, fit(i).y_pos];
        dist_vec = dist2(spot_loc,cn).^0.5;
        fit(i).dist_2_periphery = min(dist_vec);
    end
         
        
    %Debug section
    if(debug)
        figure(10);
        hold on
        mip_tmp = max(search_stack,[],3);
        imagesc(mip_tmp)
        POS = fit(i).pos - [y(1)-1,x(1)-1,z(1)-1];
        plot(POS(2),POS(1),'*r')
        hold off
        set(gca,'ydir','reverse')
        colormap gray

        figure(8);
        imshow3D(search_stack)
        hold on
        plot(POS(2),POS(1),'*r')
        hold off

        pause
    end
end 

end
