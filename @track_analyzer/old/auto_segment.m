%Nucleus segmentation automated for track_analyzer

function auto_segment( obj, params )


%Get exp info
exp_info = obj.exp_info;

%Get list of image files
img_files = dir([exp_info.nuc_seg_dir,'*.tif'])


%Define parameters for segmentation 
    seg_params = exp_info.seg_params;
    gaussSD          = seg_params(1);
    hForMaxima       = seg_params(2);
    hForMinima       = seg_params(3);
    minSizeOfNucleus = seg_params(4);
    matlabOrDipimage = seg_params(5);
    weakBorderThr1   = seg_params(6);
    weakBorderThr2   = seg_params(7);
    
%Now loop over these files and segment
parfor i = 1:length(img_files)
    

    %Test segmentation with image1
    img = imread([exp_info.nuc_seg_dir,img_files(i).name]);
    %Flatten to one stack if multiple planes
    if(size(img,3)>1)
        img = img(:,:,1);
    end
    
    %Nuclear segmentation (no plotting, non-interactive)
    [nuclearMask,~] = wahlbynucleus(img,gaussSD,hForMaxima,hForMinima,minSizeOfNucleus,matlabOrDipimage,weakBorderThr1,weakBorderThr2,0,0);
    disp(['Segmented frame: ',num2str(i,2),'/',num2str(length(img_files))])
    
    %Select areas that are within the specified range in parameters
    BW = bwareafilt( logical( nuclearMask ), params.nuc_size_range );
    

    %Save the outlines as contours
    B = bwboundaries(BW,'noholes');

    %Create data for tracking / computing:
    [refined_cells, centroids] = reorganize( B);
        
    %Create Frame obj
    frame_obj = createFrameObj( refined_cells, centroids, BW);
    
    %Now save this data
    [a,b,c] = fileparts([exp_info.nuc_seg_dir,img_files(i).name]);
    frame_file = [a,'/',b,'.mat'];
    parsave(frame_file,frame_obj)
end



%newFigure(9)
%img = imread('/home/jan/Dropbox/code_bank/matlab/cell_tracking/wahlby_testimage.tif');
%wahlbynucleus(img,1,5,3,100,1,25,30,0)



end

%Save data outside of parallel loop
function parsave(fname, frame_obj)
save(fname, 'frame_obj')
end

%Reorganize data
function [refined_cells, centroids] = reorganize(B);
for j = 1:length(B)
    %Cell contour (flip columns 1+2)
    refined_cells{j} = B{j}(:,[2,1]);

    %Calculate centroid of each object
    centroids{j} = mean(refined_cells{j});

end

end

%Generate frame structure
function frame_obj = createFrameObj( refined_cells, centroids, BW);

frame_obj.refined_cells = refined_cells;
frame_obj.centroids = centroids;
frame_obj.BW = BW;



end