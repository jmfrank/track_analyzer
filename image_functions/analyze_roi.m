%% ROI analysis function.

function data = analyze_roi(img, clean_mask, int_mask, ID, z_planes , params )

%% Nucleus. 

%First check if clean_mask and int_mask are equal.
if isequal(clean_mask,int_mask)
    
    %Erode mask a by in-boundary pixels.First pad with zeros. << forgot why
    %padding?
    mask_pd = padarray(clean_mask,[1,1],0);
    
    %If we need buffer away from segmentation edge. 
    if( params.in_buffer>0)
        se = strel('disk',params.in_buffer,8);
        mask_outter = imerode(mask_pd,se);
    else
        mask_outter = mask_pd;
    end
    
    se = strel('disk',params.in_boundary,8);
    mask_inner = imerode(mask_outter,se);
    %unpad.
    mask_inner = mask_inner(2:end-1,2:end-1);
    mask_outter = mask_outter(2:end-1,2:end-1);
    %Use some mask differences to query the ROIs. 
    roi_inner = logical( mask_outter - mask_inner);
else
    
    roi_inner = int_mask;
    
    
end

%% Cytoplasm. 
    %Now get cyto mask.
    
    %Buffer as needed
    if(params.out_buffer > 0)
        se = strel('disk',params.out_buffer,8);
        mask_inner = imdilate(clean_mask,se);
    else
        mask_inner = clean_mask;
    end
    se = strel('disk',params.out_boundary,8);
    mask_outer = imdilate(mask_inner,se);

    roi_outer = logical( mask_outer - mask_inner);
    

    %% Calculate the mean vals of nuclear and cytoplasmic regions

    %Make all masks the correct depth in 3D.
    roi_inner = repmat(roi_inner,[1,1,length(z_planes)]);
    roi_outer   = repmat(roi_outer,[1,1,length(z_planes)]);

    %Get mean/median nuclear and cytoplasm fluorescence 
    data.nuc_mean   = mean( img( roi_inner ),'omitnan' );
    data.nuc_med    = median( img( roi_inner ),'omitnan' );
    data.cyto_mean  = mean( img( roi_outer   ),'omitnan' );   
    data.cyto_med   = median( img(roi_outer),'omitnan');
    data.cell_id    = ID;
    data.local_rho  = NaN;
    data.nuc_area   = sum( clean_mask(:) );
    if(isnan(data.nuc_mean) || isnan(data.cyto_mean) || isnan(data.cyto_med) || isnan(data.nuc_med))
        disp('warning: nan found')
        disp(['cell: ',num2str(ID)]);
    end


%% DEBUGGING code for plotting contours. 
debug=0;
if debug 
    c = contourc(double(max(roi_inner,[],3)),[0.5,0.5]);
    C = C2xyz(c);
    newFigure(18)
    imshow3D(img);    
        hold on

    for i = 1:length(C)
        plot(C{i}(:,1),C{i}(:,2),'linewidth',2)
    end
end



end