%% Iterative thresholding on MCP-mNeon cells with variable expression levels. 

%I is smoothed image
%J is segmentation image. 
%real_bg is the intensity of background for image J. 
%thresh start is the starting threshold level. 
%iterations is 

function BW = iterative_thresholding(I, J, real_bg, params )


%For now, just hard code the peak range 
params.peak_range=[5:85];

%First find local maxima of I. 
BWmax = imextendedmax(I,params.Hdepth);

%Local max region props. 
stats_max = regionprops(BWmax,J,'centroid','MeanIntensity','PixelIdxList');

%Find blobs using lower bound estimation of threshold. 
BW = J >= params.thresh_start;
%Fill holes. 
BW = imfill(BW,'holes');
%Close. 
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

for i = 1:length(matched_blobs)
    %Get J values. 
    px = stats_J( matched_blobs(i) ).PixelIdxList;
    
    %Get xy range. 
    [Y,X] = ind2sub(size(BW),px);
    yrange=min(Y):max(Y);
    xrange=min(X):max(X);
    
    %Sub image of this blob. 
    sub_img = J(yrange,xrange);
    
    %Iterator until histogram is good. 
    [subBW] = iterator( sub_img, params );
    
    %Add to BW.
    BW(yrange,xrange) = BW(yrange,xrange) | subBW;
    
      
end

end


%Iteration sub-function. 
function [subBW] = iterator( img, params )

%Percentile range. 
p_range=[1:0.2:100];

%Initialize. 
thresh = params.thresh_start;
go=1;

while go
    
    %Perform analysis. 
    px = img>=thresh;
    vals = norm_data(img(px));
    P =  smooth(prctile(vals,p_range),12);
    dP = diff(P);
    [~,LOCS] = findpeaks(dP,p_range(2:end),'MinPeakProminence',0.003);
    these_pks = find(LOCS > params.peak_range(1) & LOCS < params.peak_range(end));

    %Check peaks. 
    if isempty(these_pks)
        
        go=0;
        
    else
        
        %Increasing threshold. 
        thresh = thresh + params.increment;
        
    end

end

subBW = img >= thresh;

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

