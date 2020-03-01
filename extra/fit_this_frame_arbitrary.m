%Pass a 3D frame and find spots. Uses the stats object pixels to query the
%local background to give blob measurements. 

function fit = fit_this_frame_arbitrary( img, stats, BW, params )

debug = 0;
%Get dimensions of stack
[l,w,h] = size(img); 

% if mask is 2D, expand to 3D
if length(size(BW)) == 2
    BW = repmat(BW,[1,1,h]);
end

% number of centroids. 
N = length(stats);

%Get vectors of peak centers and mean-intensities
centers = cat(1,stats.Centroid);

%Empty structure for fitting measurements. Input is the number of
%centroids. Will need to remove empty entries later on. 
fit = gen_fit_struct( N );

% We will use image dilation to sample the local background. Need to ignore
% pixels that are not part of cell.
bounds=params.bg_mask;
B = strel('sphere',bounds(2));
A = strel('sphere',bounds(1));
bad_list=[];
for i = 1:N
        
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
    img_sub_stack = double(img(y,x,z));
    
    % Substack of this centroid's mask.
    img_thresh = false(size(BW));
    img_thresh( stats(i).PixelIdxList ) = 1;
    fg_sub_stack = img_thresh(y,x,z);
    
    % Get roi in cell BW.
    bg_sub_stack = BW(y,x,z);
    
    % generate the background mask
    this_bg = imdilate(fg_sub_stack, B) - imdilate(fg_sub_stack, A);
    
    % Now multiply by cell mask.
    this_bg = logical(this_bg.*bg_sub_stack);
    
    % Calculate mean bg intensityies.
    bg_int = img_sub_stack(this_bg);
    mean_bg = mean( bg_int );
    bg_std  = std(bg_int);
    
    % Create threshold for foreground pixels based on bg+ std*X. This
    % allows more accurate measurement of blob wrt local background
    % intensities. 
    %fg_thresh = mean_bg + params.STD_threshold*bg_std;
    %fg_mask = img_sub_stack >= fg_thresh & imdilate(fg_sub_stack, A).*bg_sub_stack;
    %S = regionprops(fg_mask,'Area','PixelIdxList');
    % Take maximum size region. 
    %[~,idx] = max([S.Area]);
    %this_fg = S( idx ).PixelIdxList;
    %this_fg_mask = false(size(fg_mask));
    %this_fg_mask(this_fg) = 1;
    this_fg = find(fg_sub_stack);
    this_fg_mask = fg_sub_stack;
    
    %Background subtracted foreground values
    fg_int = img_sub_stack( this_fg ) - mean_bg;
    
    %foreground coordinates
    [Y,X,Z] = ind2sub(size(img_sub_stack),this_fg);
    
    %Center of mass. 
    M = sum(fg_int);
    R = sum([Y,X,Z].*fg_int,1)./M;
    
    %Fill up fit variable. 
    fit(i).pos      = R + [y(1)-1,x(1)-1,z(1)-1]; %(Shift back to position)
    fit(i).mean_int = mean(fg_int);    
    fit(i).real_bg  = mean_bg;
    fit(i).sum_int  = sum(fg_int);
    fit(i).snr      = max(fg_int)/std(bg_int);
    
    if isempty(fit(i).snr) | isinf(fit(i).snr) | isnan(fit(i).snr)
        
        %fit(i).sum_int = [];
        bad_list = [bad_list, i];
    end
    %Assign fit to cell
    fit(i).cell_id  = stats(i).assignment;
    
    % Size of blob
    fit(i).size     = length(this_fg);
    
    % Pixel idx of blobs. 
    pixels = false(size(img));
    pixels(y,x,z) = this_fg_mask;
    fit(i).PixelIdxList = find(pixels);
    
    %Debug section
    if(debug)
        POS = fit(i).pos - [y(1)-1,x(1)-1,z(1)-1];
        figure(8);
        imshow3D(img_sub_stack)
        hold on
        plot(POS(2),POS(1),'*r')
        hold off
    end
    
end 


% Remove bad fits. 
keep = 1:length(fit);
keep = setdiff(keep,bad_list);
fit = fit(keep);


end

function empty_structure = gen_fit_struct( N )

%Fields
empty_structure(N) = struct('pos',[],'mean_int',[],...
    'real_bg',[],'sum_int',[],'cell_id',[],'size',[],'PixelIdxList',[]);

end



