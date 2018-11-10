%Open an image. Max-proj. Display time lapse as '3d' image. Ask user for
%ROI. Get ROI intensity trace. Useful for getting mean-background signal. 


function obj = go_get_roi(obj)


debug = 0;

%Initiate track. 
frap_track = [];

%Load whole image. 
bf_img = bfopen( obj.exp_info.img_file );
reader = bfGetReader(obj.exp_info.img_file);

%Now run through the rest of the time points and track. 
T = reader.getSizeT;

disp('Loading image...')
for t = 1:T
    stack = get_stack(reader,t,1);
    %Append to IMG. 
    IMG(:,:,t)=max(stack,[],3);

end
reader.close();

%Imshow3D. 
figure(2);
imshow3D(IMG);

%Ask for ROI. 
h = impoly('Closed',1);

%Get positions of polygon. 
vec = h.getPosition;

%Convert to mask. 
MASK = poly2mask(vec(:,1),vec(:,2),size(IMG,1),size(IMG,2));

%Get pixel idx. 
PixelIdx = find(MASK);

%Now make a FRAP object. 
frap = frap_analyzer( obj.exp_info );
roi.PixelIdx = PixelIdx;
roi.type     = 'polygon';
roi.img_type = 'img';

%Calculate mean trace with open image. 
for i = 1:T
    this_img = IMG(:,:,i);
    int_trace(i) = mean( this_img( PixelIdx ) );
end
roi.trace = int_trace;
frap = frap.addROI( roi );

%Now add frap object to obj. 
obj.frap = frap;

disp('finished.')
close(gcf);
end



%Get fits based on filtering detection
function fit = get_fits(stack, params)

mip = max(stack,[],3);
h = fspecial('log',[15,15],8); %,[5,5],params.log_sigma(1));
img_filter = -imfilter(mip,h,'symmetric', 'conv');

img_bw = img_filter >= params.thresh;
stats = regionprops(img_bw,'Centroid','Area');
A = num2cell( ones(length(stats),1));
[stats.assignment] = A{:};    
if(isempty(stats))
    fit = [];
else
    fit = fit_this_frame_centroid_2D( mip, stats, params );
end

end