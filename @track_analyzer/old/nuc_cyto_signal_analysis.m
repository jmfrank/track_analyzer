% Calculate mean yap signal (could be arbitrary signal based on image).
% Uses the roi2poly tool. Maybe faster?
%Now computes cyto to nuclear ratio by taking a mean of a thin ring outside
%of the contour

%***Updated to use structures for calculations. 

%***2-27-18 Now we will search through the stack to find the best plane(s) for
%calculating N/C signals. We will find which plane gives the maximum
%nuclear intensity for the pre-calculated mask. Then take that plane +/- a
%few other planes to do MIP and nuc/cyto calculations

function obj = nuc_cyto_signal_analysis(obj, params)


debug = 0;
if(debug)
    figure(7)
    clf(7)
    %set(gcf,'color','w')    
    color_vec  = hsv(length(obj.tracks));
    color_vec  = color_vec(randperm(size(color_vec,1)),:);
end

Z = obj.exp_info.z_planes;
T = obj.exp_info.t_frames;


%UPDATED: USING BIOFORMATS READER
%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = bfGetReader(obj.exp_info.img_file);
series = 1;

%Get the image size of this series. 
size_x = reader.getSizeX;
size_y = reader.getSizeY;


%Generate pixel location matrices
xgv=[1:size_x];
ygv=[1:size_y];
[X,Y]=meshgrid(xgv,ygv);

%Get segmentation file names. 
seg_files = obj.get_frame_files;

for t = 1:T
    display(['Analyzing mean signal of frame: ',num2str(t)])
    tic
    %Create empty image to fill up
    msk_img = zeros(size_y,size_x,Z);
    sig_img = zeros(size_y,size_x,Z);
    
    %Get the bio-formats image index corresponding to this z-stack:
    for i = 1:Z
        this_plane = reader.getIndex(i-1,params.seg_channel-1,t-1)+1;
        msk_img(:,:,i) = bfGetPlane(reader,this_plane);
        
        this_plane = reader.getIndex(i-1,params.signal_channel-1,t-1)+1;
        sig_img(:,:,i) = bfGetPlane(reader,this_plane);
    end
    
        
    %Load the Frame_obj file associated with this image
    load(seg_files{t});
    
    %Loop over detected cells
    try
        n_cells = length(frame_obj.refined_cells);
        seg_type = 'refined';
    catch
        n_cells = length(frame_obj.PixelIdxList);
        seg_type = 'Px';
    end
    
    %Create empty structure for nuclear / cytoplasm calculations
    data = gen_data_struct( n_cells );
        
    
    
    for j = 1:n_cells
        
        %Get a 2D mask of this cell, depending on segmentation format....
        switch seg_type
            case 'refined'
                %Contour
                x     = frame_obj.refined_cells{j}(:,1);
                y     = frame_obj.refined_cells{j}(:,2);
                %Get mask of this cell
                roi_bw_nuc = poly2mask(x,y,size(sig_img,1),size(sig_img,2));  
            case 'Px' 
                BW = zeros(size(sig_img,1),size(sig_img,2));
                px = frame_obj.PixelIdxList{j};
                BW(px) = 1;
                %Create contour. 
                C = contourc(BW,[0.5,0.5]);
                [X,Y] = C2xyz(C);
                %Sometimes, there multiple contours. Loook for largest one?
                if(length(X)>1)
                    S  = cellfun('length',X);
                    [~,id] = max(S);
                    x = X{id};
                    y = Y{id};
                else                    
                    x = X{1};
                    y = Y{1};
                end
                %Create mask. 
                roi_bw_nuc = logical(BW);
        end
        
        %% Figure out optimal z-section. 
        %loop over planes
        for z = 1:Z
            this_plane = msk_img(:,:,z);
            mean_z_val(z) = mean( this_plane(roi_bw_nuc) );
        end
        %Find index of max value <-i.e. the 'best' plane
        [~,idx] = max(mean_z_val);
        %Make a range of planes to use based on params.z_range
        z_planes = idx-params.z_range: idx+params.z_range;
        %Need to adjust z-planes if they are outside of img
        if(z_planes(end) > Z)
            sel = z_planes <= Z;
            z_planes = z_planes(sel);
        end
        if( z_planes(1) < 1 )
            sel = z_planes >= 1;
            z_planes = z_planes(sel);
        end

       %% Creating ROIs for nuclear and cytoplasm 
            %Mean
            mean_x = mean(x);
            mean_y = mean(y);

            %Make outer and inner contours. 
            vec = 
            
            %Make outer contour
            clear x_out y_out y_in x_in
            for k = 1:length(x)
                vec = [x(k)-mean_x,y(k)-mean_y];
                vec_norm = vec ./ sqrt( sum(vec.^2) );
                x_out(k) = x(k) + vec_norm(1)*params.out_boundary;
                y_out(k) = y(k) + vec_norm(2)*params.out_boundary;
                x_in(k)  = x(k) - vec_norm(1)*params.in_boundary;
                y_in(k)  = y(k) - vec_norm(2)*params.in_boundary;
            end


            %Determine if each point is within contour
            if(step.gpu)
                roi_bw_nuc = poly2mask_gpu(x_in,y_in,size(sig_img,1),size(sig_img,2));
                %Outer contour
                roi_bw_in  = poly2mask_gpu(x,y,size(sig_img,1),size(sig_img,2));
                roi_bw_out = poly2mask_gpu(x_out,y_out,size(sig_img,1),size(sig_img,2));
            else
                roi_bw_nuc = poly2mask(x_in,y_in,size(sig_img,1),size(sig_img,2));
                %Outer contour
                roi_bw_in  = poly2mask(x,y,size(sig_img,1),size(sig_img,2));
                roi_bw_out = poly2mask(x_out,y_out,size(sig_img,1),size(sig_img,2));
            end
            %Difference contour
            roi_diff = logical(roi_bw_out - roi_bw_in);
        %% Calculate the mean vals of nuclear and cytoplasmic regions
            
            %Get img based on z_range
            img_dbl = double(sig_img(:,:,z_planes));
            
            %Make all masks the correct size in 3D.
            %roi_bw_nuc = repmat(roi_bw_nuc,[1,1,length(z_planes)]);
            %roi_diff   = repmat(roi_diff,[1,1,length(z_planes)]);
            
            %Get mean nuclear and cytoplasm fluorescence 
            data(j).nuc_mean   = mean( img_dbl( roi_bw_nuc ) );
            data(j).nuc_med    = median( img_dbl( roi_bw_nuc ) );
            data(j).cyto_mean  = mean( img_dbl( roi_diff   ) );   
            data(j).cyto_med   = median( img_dbl(roi_diff));
            data(j).cell_id    = j;
            
            if(isnan(data(j).nuc_mean))
                
                pause
            end
            
                    
        %% Plot during debugging
        if(debug)
            figure(12)
            subplot(1,3,1)
            imshow(roi_bw_nuc);
            subplot(1,3,2)
            mip_tmp = max( img_dbl,[],3);
            imagesc(mip_tmp);
            subplot(1,3,3)
            imagesc(mip_tmp);
            hold on
            colormap gray
            plot(x_in,y_in,'color',color_vec(j,:),'linewidth',1)
            plot(x_out,y_out,'color','y','linewidth',1)
            hold off
            pause
        end
    end
    
    %Now save the frame_obj
    frame_obj.data = data;
    save(seg_files{t},'frame_obj','-append')

    if(debug);pause;if(i < length(obj.img_files));clf(7);end;end
    toc
end

       
if(debug);hold off;end;

obj.nuc_cyto_calc = true;

%Now 
end


%% Generate an empty structure for keeping track of nuclear and cytoplasmic signals
function data = gen_data_struct( N )

%Fields
data(N) = struct('nuc_mean',[],'cyto_mean',[],'local_rho',[],'cell_id',[],'nuc_med',[],'cyto_med',[]);

end