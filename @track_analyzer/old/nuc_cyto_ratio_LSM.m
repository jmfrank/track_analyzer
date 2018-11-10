% Calculate mean yap signal (could be arbitrary signal based on image).
% Uses the roi2poly tool. Maybe faster?
%Now computes cyto to nuclear ratio by taking a mean of a thin ring outside
%of the contour

%***Updated to use structures for calculations. 

%***2-27-18 Now we will search through the stack to find the best plane(s) for
%calculating N/C signals. We will find which plane gives the maximum
%nuclear intensity for the pre-calculated mask. Then take that plane +/- a
%few other planes to do MIP and nuc/cyto calculations



function obj = nuc_cyto_ratio_LSM(obj, params)


debug = 0;
if(debug)
    figure(7)
    clf(7)
    %set(gcf,'color','w')    
    color_vec  = hsv(length(obj.tracks));
    color_vec  = color_vec(randperm(size(color_vec,1)),:);
end

Z = obj.exp_info.z_frames;
T = obj.exp_info.t_frames;

%For quickly querying the signal location within contour
img = tiffread(obj.exp_info.img_file,1);
%Generate pixel location matrices
[img_y,img_x] = size(img.data{1});
xgv=[1:size(img.data{1},2)];
ygv=[1:size(img.data{1},1)];
[X,Y]=meshgrid(xgv,ygv);


for t = 1:T
    display(['Analyzing mean signal of frame: ',num2str(t)])
    %load signal image
    
    switch f_type
        case '.lsm'
            
            %Calculate indices. Skip every other b/c thumbnails? 
            start = Z.*(t-1)*2 + 1;
            indices  = start:2:start + 2*Z-1;
            %Grab this Z-stack
            IMG = tiffread( obj.exp_info.img_file, indices );

            %Now project
            sig_img = zeros(img_y,img_x,Z);
            for z = 1:Z
                %Signal image
                sig_img(:,:,z) = IMG(z).data{params.signal_channel};
                %Mask image
                msk_img(:,:,z) = IMG(z).data{params.seg_channel};

            end
            
        case '.tif'
            %UPDATE THIS SECTION
            
    end
    
    
        
    %Load the Frame_obj file associated with this image
    fname = ['frame_',sprintf('%04d',t),'.mat'];
    
    load([obj.exp_info.nuc_seg_dir,fname],'frame_obj');
    
    %Create empty structure for nuclear / cytoplasm calculations
    data = gen_data_struct( length( frame_obj.refined_cells) );
    
    %Loop over detected cells
    for j = 1:length(frame_obj.refined_cells)
        
        %Contour
        x     = frame_obj.refined_cells{j}(:,1);
        y     = frame_obj.refined_cells{j}(:,2);

        %% Now find plane that gives maximum of the DNA channel
            %Get mask of this cell
            roi_bw_nuc = poly2mask(x,y,size(sig_img,1),size(sig_img,2));    
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
            roi_bw_nuc = poly2mask(x_in,y_in,size(sig_img,1),size(sig_img,2));
            %Outer contour
            roi_bw_in  = poly2mask(x,y,size(sig_img,1),size(sig_img,2));
            roi_bw_out = poly2mask(x_out,y_out,size(sig_img,1),size(sig_img,2));
            %Difference contour
            roi_diff = logical(roi_bw_out - roi_bw_in);
        %% Calculate the mean vals of nuclear and cytoplasmic regions
            
            %Get img based on z_range
            img_dbl = double(sig_img(:,:,z_planes));
            
            
            %Get mean nuclear and cytoplasm fluorescence 
            data(j).nuc_mean   = mean( img_dbl( roi_bw_nuc ) );
            data(j).cyto_mean  = mean( img_dbl( roi_diff   ) );   
            data(j).cell_id    = j;
        
            nuc_cyto_ratio=data(j).nuc_mean/data(j).cyto_mean;
            
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
    save([obj.exp_info.nuc_seg_dir,fname],'frame_obj','-append')

    if(debug);pause;if(i < length(obj.img_files));clf(7);end;end
end

       
if(debug);hold off;end;

obj.nuc_cyto_calc = true;

end


%% Generate an empty structure for keeping track of nuclear and cytoplasmic signals
function data = gen_data_struct( N )

%Fields
data(N) = struct('nuc_mean',[],'cyto_mean',[],'local_rho',[],'cell_id',[]);


end