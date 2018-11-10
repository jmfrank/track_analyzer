% Uses the masks in the obj.mask field to calculate nuclear / cyto ratio.
% Will overwrite the tracks and nuc_cyto_data fields so that they match the
% masks. 

%Updated to use structures for calculations. 
function obj = nuc_cyto_ratio_masks(obj, params, step)


debug=0;

if nargin <2
    params=struct;
end
params = default_params(params);

%Number of cells is based on the number of mask sets
n_cells = length(obj.masks);
%Make n_cells empty tracks
obj.tracks = cell(n_cells,1);
%Clear out data. 
obj.nuc_cyto_data = struct;
%Loop over time and fill out tracks and nuc_cyto_data
data_count = 0;
    
%Create empty structure for nuclear / cytoplasm calculations. 
data = gen_data_struct( n_cells );


%Look too see if a max-p image exists. 
if( exist( obj.exp_info.max_p_img,'file' ) )

    %Just one image. 
    fname = obj.exp_info.max_p_img;

    [reader,X,Y,Z,C,T] = bfGetReader(fname);
    this_img = bfopen(fname);

    frame_range = [1:T]';
    img_val = 1.*ones(length(frame_range),1);
    %Add to frame2img. 
    frame2img=[frame_range,frame_range, img_val];

    Img = zeros(Y,X,T);
    for t = 1:T
        plane = get_planesZCT(reader,Z,params.channel,t);
        Img(:,:,t) = this_img{1}{plane,1};
    end
    
    %Close reader
    reader.close();
    
else
    error('missing max-p img');

end



for t = 1:T
    
    display(['Analyzing mean signal of frame: ',num2str(t)])

    %This time-point 
    img_dbl = double(Img(:,:,t));
    
    if(debug);imagesc(img);colormap('gray');end;

    %Loop over detected cells
    for j = 1:n_cells
        %update counter
        data_count = data_count + 1;
        
        %Get masks
        roi_nuc = obj.masks(j).nuc_mask;
        roi_cyto = obj.masks(j).cyto_mask;
        
        %Get mean nuclear and cytoplasm fluorescence 
        data(data_count).nuc_mean   = mean( img_dbl(  roi_nuc    ) );
        data(data_count).cyto_mean  = mean( img_dbl(  roi_cyto   ) );   
        data(data_count).cell_id    = j;
        
        %get centroid of nuclear blob
        stats = regionprops(roi_nuc,'centroid');
        x = stats(1).Centroid(1);
        y = stats(1).Centroid(2);
        %Make track entry
        obj.tracks{j}(t,:) = [t,NaN,x,y, data_count ];
        
        %Plot during debugging
        if(debug)
            figure(12)
            subplot(1,3,1)
            imshow(roi_bw_nuc);
            subplot(1,3,2)
            imshow(roi_diff);
            subplot(1,3,3)
            imagesc(img);
            colormap gray
            hold on
            plot(x_in,y_in,'color',color_vec(j,:),'linewidth',1)
            plot(x_out,y_out,'color','y','linewidth',1)
            hold off
            pause
        end
    end
    
    if(debug);pause;if(t < length(obj.img_files));clf(7);end;end
end

%Figure out channel name. 
channel_str = ['channel_',pad( num2str(params.channel),2,'left','0')];
obj.nuc_cyto_data.(channel_str) = data;

       
if(debug);hold off;end;

obj.nuc_cyto_calc = true;
obj.tags = []; %Clear tags


%Look for pre-image data. 
if step.analyze_pre_image

    %Open up the pre-img. 
    IMG = bfopen(obj.exp_info.pre_img_file);
    [reader,X,Y,Z,C,T] = bfGetReader(obj.exp_info.pre_img_file);

    planes = get_planesZCT(reader,Z,params.pre_img_channel,1);
    stack = zeros(Y,X,Z);

    for z =1:Z
        stack(:,:,z) = IMG{1}{planes(z),1};
    end

    img = max(stack,[],3);

    %Check if there's a frame_obj. 
    %Save frame_obj. 
    pre_dir = obj.exp_info.pre_img_nuc_seg_dir;
    pre_img_frame=[pre_dir,'/frame_0001.mat'];
    
    if exist(pre_img_frame,'file')
        load(pre_img_frame);
    else
        %Check there are enough cells in frame_obj. 
        n_cells = length(obj.masks);

        %Get roi info. 
        exp_info.foreground_file = obj.exp_info.img_file;
        exp_info.background_file=[];
        frap = frap_analyzer(exp_info);
        frap = frap.getROIinfo;

        %Close reader
        reader.close();
        
        if step.use_roi
            %loop over roi's. 
            for r = 1:length(frap.ROI)
                roi=frap.ROI(r);
                BW = poly2mask(roi.x_pos,roi.y_pos,Y,X);
                frame_obj.PixelIdxList{r} = find(BW);
                [y,x] = ind2sub( size(img), find(BW));
                frame_obj.centroid{r}  = [mean(x),mean(y)];
            end
            
        else
            clf;
            figure(7);
            clf;
            imshow3D(img);
            hold on
            for  r = 1:length(frap.ROI)
                roi= frap.ROI(r);
                plot(roi.x_pos,roi.y_pos,'--','linewidth',2,'color',[1,0.5,0.5,0.5]);
            end
            hold off
            %Loop over cells, collect a user-generated ROI. 
            for i = 1:n_cells       
                %First create nuclear ROI mask;
                this_mask = roipoly;

                disp(['captured roi # ',num2str(i)])
                frame_obj.PixelIdxList{i} = find(this_mask);
                [y,x] = ind2sub( size(img), find(this_mask));
                frame_obj.centroid{i}  = [mean(x),mean(y)];
            end
        end
        
        save(pre_img_frame,'frame_obj');

    end

    %Error if different number of masks and objects. 
    if length(obj.masks) ~= length(frame_obj.centroid)
        error('missing a mask or centroid')
    end
    
    %Match masks to hand-drawn cells. 
    pre_centroids = cell2mat(frame_obj.centroid');
    for m = 1:length(obj.masks)
        roi = regionprops( obj.masks(m).nuc_mask,'centroid');
        mask_centroids(m,:) = roi.Centroid;
    end
    D = pdist2(pre_centroids,mask_centroids);
    [~,idx] = min(D);
    frame_obj.PixelIdxList = frame_obj.PixelIdxList(idx);
    frame_obj.centroid = frame_obj.centroid(idx);
    

    
    
    %Loop over cells in pre-img and measure intensity properties. 
    for i = 1:length(frame_obj.centroid)
        %Now query the ROI for the pre-image. 
        pre_data(i).nuc_mean = mean( img( frame_obj.PixelIdxList{i}));
        pre_data(i).nuc_sum  = sum( img( frame_obj.PixelIdxList{i}));
        pre_data(i).nuc_med = median( img( frame_obj.PixelIdxList{i}));
    end 
   
    %Append the data to a pre-img field of nuc_cyto_data. 
    obj.nuc_cyto_data.pre_img = pre_data;
end

    
    
end


%% Generate an empty structure for keeping track of nuclear and cytoplasmic signals
function data = gen_data_struct( N )

%Fields
data(N) = struct('nuc_mean',[],'cyto_mean',[],'local_rho',[],'cell_id',[]);


end