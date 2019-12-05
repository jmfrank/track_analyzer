%Pass a 3D frame and find spots. Uses 'region_props' data to make fits. 

function fit = fit_this_frame_centroid( img, stats, params )

debug = 0;
%Get dimensions of stack
[l,w,h] = size(img); 

%Get vectors of peak centers and mean-intensities
centers = cat(1,stats.Centroid);

%Empty structure for fitting measurements. Input is the number of
%centroids. Will need to remove empty entries later on. 
fit = gen_fit_struct( length(stats) );

%Make a special shell-mask for querying the local background. 
    bounds=params.bg_mask;
    B = strel('sphere',bounds(2));
    A = strel('sphere',bounds(1));
    
    fgMASK = padarray(A.Neighborhood,(params.search_radius-bounds(1))*[1,1,1]);
    bgMASK = padarray(B.Neighborhood,(params.search_radius-bounds(2))*[1,1,1])...
        & ~fgMASK;

for i = 1:size(centers,1)
        
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
    this_fg = fgMASK(sel_y,sel_x,sel_z);
    this_bg = bgMASK(sel_y,sel_x,sel_z);
    
    %Background calculation. 
    bg_int = search_stack( this_bg );
    mean_bg = mean(bg_int);
    
    %Foreground index for search stack
    fg_ind = find(this_fg);
    
    %Background subtracted foreground values
    fg_int = search_stack( fg_ind ) - mean_bg;
    
    %foreground coordinates
    [Y,X,Z] = ind2sub(size(search_stack),fg_ind);
    
    %Center of mass. 
    M = sum(fg_int);
    R = sum([Y,X,Z].*fg_int,1)./M;
    
    %Fill up fit variable. 
    fit(i).pos      = R + [y(1)-1,x(1)-1,z(1)-1];
    fit(i).mean_int = mean(fg_int);    
    fit(i).real_bg  = mean_bg;
    fit(i).sum_int  = sum(fg_int);
    fit(i).snr      = max(fg_int)/std(bg_int);
    %Assign fit to cell
    fit(i).cell_id  = stats(i).assignment;
           
        
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

function empty_structure = gen_fit_struct( N )

%Fields
empty_structure(N) = struct('pos',[],'mean_int',[],...
    'real_bg',[],'sum_int',[],'cell_id',[]);

end



