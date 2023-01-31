%% Iterative thresholding on MCP-mNeon cells with variable expression levels. 

%I is smoothed image
% Perform local thresholding. This version first looks for local max. Then
% we generate a local threshold from these local maxima. Just use 0.5x the
% local max intensity seems to work well. We also need a minimum intensity
% to ignore local maximums in background. 

function BW = local_thresholding(J, params )


%For now, just hard code the peak range 
params.peak_range=[5:88];

%First find local maxima of I. 
BWmax = imextendedmax(J,params.Hdepth);

%Local max region props. 
stats_max = regionprops(BWmax,J,'centroid','MeanIntensity','PixelIdxList');

% Use the mean intensity of blobs to create an estimate of the best local
% intensity threshold. 
C = cat(1,stats_max.Centroid);
ints = cat(1,stats_max.MeanIntensity);
img_size = size(J);
[xq,yq] = meshgrid(1:img_size(1),1:img_size(2));
F = scatteredInterpolant( C(:,1),C(:,2),ints, 'natural','nearest');
vq = F(xq,yq);

% Generate first estimate from vq and starting threshold.
BW = J >= 0.5.*vq;


%Find blobs using lower bound estimation of threshold. 
BW = J >= params.thresh_start;
%Fill holes. 
BW = imfill(BW,'holes');
%Close. BWmax
stats_J = regionprops(BW,'centroid','PixelIdxList','Area');
%Filter out J blobs that are too small. 
sel = [stats_J.Area] >= 0.8*params.AbsMinVol;
stats_J=stats_J(sel);


%Match the regional maxima to blobs in BW. 
ctr_maxima = cat(1,stats_max.Centroid);
ctr_J      = cat(1,stats_J.Centroid);
D = pdist2(ctr_maxima,ctr_J);
Dsel = D <= 3*sqrt(params.AbsMaxVol/pi);

%Loop over regional max. 
matched=zeros(1,length(stats_max));
for i = 1:length(stats_max)
    
    %Are there any blobs nearby. 
    blob_ids = find( Dsel(i,:) );
    
    %Skip if empty. 
    if isempty( blob_ids )
        continue
    end
    
    %Loop over blob_ids. 
    for b = blob_ids
        
        %Query if blob encompasses most (or all?) of regional max. 
        unmatched_pixels = setdiff( stats_max(i).PixelIdxList, stats_J( b ).PixelIdxList );
        
        if length(unmatched_pixels) <= 0.09*length(stats_max(i).PixelIdxList)
            matched(i) = b;
        end
    end
end

%Unique blobs. 
matched_blobs = unique( matched( matched > 0 ) );
unmatched_peaks = find(~matched);
BWclean = false(size(BW));
px = cat(1,stats_J( matched_blobs ).PixelIdxList );
BWclean(px) = 1;

%Loop over matched blobs. Evalulate histogram, increase threshold until
%critereum met. 
BW=false(size(BW));
%Dilating structural element. Square is faster than disk. 
params.dilator_size=5;
dilator = strel('square',params.dilator_size);

%%% Note: somehow this happened to work for 3D input image. Even though
%%% the syntax looks like 2D (i.e. only specifying Y,X range). Could be an
%%% issue with more modifications. 
for i = 1:length(matched_blobs)
    %Get J values. 
    px = stats_J( matched_blobs(i) ).PixelIdxList;
    
    %Build mask of this cell. 
    this_cell_mask= false(size(BW));
    this_cell_mask(px) = 1;
    %Dilate a bit. 
    this_cell_mask = imdilate(this_cell_mask,dilator);
    %Get range.
    [Y,X] = find(this_cell_mask);
    yrange=unique(Y);
    xrange=unique(X);    
    
    %mask. 
    mask = false(size(BW));
    mask(yrange,xrange) = 1;
    %Find pixels part of other matched blobs. 
    sel = matched_blobs ~= matched_blobs(i);
    not_these_px = cat(1,stats_J(matched_blobs(sel)).PixelIdxList);
    mask(not_these_px) = 0;
    
    % Change
    % Sub image of this blob. 
    int_vals = J(mask);
    
    %Iterator until histogram is good. 
    thresh = iterator( int_vals, params );
    
    %Apply threshold to this mask. 
	subImg = J(yrange,xrange);
    subBW =  subImg >= thresh & mask(yrange,xrange);
        
    %Add to BW.
    BW(yrange,xrange) = BW(yrange,xrange) | subBW;
    
end

end


%Iteration sub-function. 
function thresh = iterator( int_vals, params )

%Percentile range. 
p_range=[1:0.2:100];

%Initialize. 
thresh = params.thresh_start;
params.prom_threshold=0.0003;
params.maxPeakP = 0.5;
go=1;

while go
    
    %Perform analysis. 
    px = int_vals>=thresh;
    vals = norm_data(int_vals(px));
    %Need to smooth the percentile before taking derivative. 
    P =  smooth(prctile(vals,p_range),12);
    %take derivative and smooth again to remove weak peaks. 
    dP = smooth(diff(P),20);
    [~,LOCS,~,PROM] = findpeaks(dP,p_range(2:end),'MinPeakProminence',params.prom_threshold);
    [~,idx_LOCS] = intersect(p_range,LOCS);
    
    %Check that the peaks satisfy histogram criterium. 
    these_pks = find(LOCS > params.peak_range(1) & LOCS < params.peak_range(end) & P(idx_LOCS)' <= params.maxPeakP );
   
    %Check peaks. 
    if isempty(these_pks)
        
        go=0;
        
    else
        
         
        LOCS = LOCS(these_pks);
        PROM = PROM(these_pks);

        %Find peak of greatest prominance.
        [~,pk_idx] = max(PROM);
        
        pct_val = LOCS(pk_idx);
        
        %Increasing threshold. 
        thresh = prctile(int_vals(px),pct_val);
        go=0;
    end
    
%Debugging.
    
%     figure(2);
%     plot(p_range,P,'linewidth',3)
%     ylabel('Intensity (normalized)')
%     yyaxis right
%     plot(p_range(1:end-1), dP,'linewidth',3)
%     ylabel('dI / dP')
%     std_figure_settings

end

%subBW = img >= thresh;

end


%% Plotting functions useful for debugging. 
function plot_p_dist( p_range,P )

newFigure(4);
plot(p_range,P,'linewidth',2);

end



%Plotting text over image
function plot_text( stats, blob_ids )

delete(findobj(gca,'Type','text'))
hold on
for i = 1:length(blob_ids)
    ctr = stats( blob_ids(i) ).Centroid;
    text(ctr(1),ctr(2),num2str( blob_ids(i)),'color','g')
end
hold off

end

%% debugging. 

% plot_text( stats_J( matched_blobs ), 1:length(matched_blobs))

%Ploting peaks. 
% figure(67); plot(p_range(2:end),dP)