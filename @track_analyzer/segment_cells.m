%Simple version of data_doctor cell segmentation. 

function obj = segment_cells(obj, params, step)

debug= step.debug;



Z = obj.exp_info.z_planes;
T = obj.exp_info.t_frames;

%UPDATED: USING BIOFORMATS READER

%open the whole image
%IMG = bfopen(obj.exp_info.img_file);

%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = bfGetReader(obj.exp_info.img_file);
series = 1;

%Get the image size of this series. 
size_x = reader.getSizeX;
size_y = reader.getSizeY;

%Loop over time and segment each frame. 
for t = 1:T
    
    %Clear out frame_obj structure
    frame_obj = struct();

    %Create empty image to fill up
    img = zeros(size_y,size_x,Z);
    
    %Get the bio-formats image index corresponding to this z-stack:
    for i = 1:Z
        this_plane = reader.getIndex(i-1,params.seg_channel-1,t-1)+1;
        img(:,:,i) =bfGetPlane(reader,this_plane);
    end
      
    %Max project
    mip = max(img,[],3);
    
    %If requested, save the mip image of the first frame to test
    %segmentation parameters on data-doctor
    if(step.MakeTmpImage & t==1)
       out_file = [obj.exp_info.nuc_seg_dir,'/frame_1.tif'];
       imwrite(uint16(mip),out_file);       
    end
 
    %Prefilter
    if(params.pre_gauss > 0)
        mip_filt = imgaussfilt_py( mip, params.pre_gauss);
    end
    
    if(params.dog(1) > 0)
        %Filter mip image
        mip_filt = imgaussfilt_py( mip_filt,params.dog(1)) - imgaussfilt_py(mip_filt, params.dog(2));
    end
    
    %% Display in debug
    if(debug)
        clf
        figure(7)
        imshow3D_filter(mip_filt);
        figure(8)
        imagesc(mip_filt)
        colormap gray
        hold on
        contour(mip_filt,[params.thresh,params.thresh],'color','w');
        set(gca,'ydir','reverse')
        colorbar
    end
    
    %Get contours
    C = contourc(mip_filt, [params.thresh,params.thresh]);
    cells = C2xyz(C);
    %Get areas 
    areas = contour_areas( cells );
    
    %Check for weird boundary effects. 
    %cells = check_for_edge_effects( cells, areas, size(mip));
    
    %Filter areas 
    sel = areas >= params.nuc_size_range(1) & areas <= params.nuc_size_range(2);
    cells = cells(sel);
    areas = areas(sel);
    
    %Remove cells that are interior to each other
    %First calculate distance matrix based on centroids
    centroids = cell2mat(cellfun(@(x) mean(x),cells,'UniformOutput',false)');
    D = pdist( centroids); 
    clean_ids = rm_duplicate_centroids(D,areas,size(centroids,1),cells, params.rm_cont_dist);

    %Create empty image
    bw = zeros( size(mip,1),size(mip,2));

    %Counter
    c = 0;
    %Loop over potential cells
    for j = clean_ids

       %temporary cell
       tmp_cell = cells{j};

       %Get mean intensity of this cell from raw image
       tmp_bw = poly2mask( cells{j}(:,1),cells{j}(:,2),size(mip,1),size(mip,2));
       tmp_int = mean( mip( tmp_bw ) ); %<< mean intensity of this cell
       %Get pixel locations of this mask
       PixelIdxList = find(tmp_bw);
       

       %Use threshold to ensure this is a fluorescence positive cell...
       if tmp_int >= params.nuc_int_thresh 
        c=c+1;
        frame_obj.refined_cells{c} = tmp_cell;
        frame_obj.centroids{c} = mean( tmp_cell );
        frame_obj.PixelIdxList{c} = PixelIdxList;
        bw = bw + tmp_bw;
       end
    end
    
    %Add the BW image to frame_obj
    frame_obj.BW = bw;
    %Add the max-int projection to frame_obj
    frame_obj.mip = mip;
    
    if(debug)
        figure(9)
        clf
        imagesc(mip);
        colormap gray
        hold on
        for c = 1:length(frame_obj.refined_cells)
            plot(frame_obj.refined_cells{c}(:,1),frame_obj.refined_cells{c}(:,2),'r')
        end
        pause 
    end
    

    
    %% Now save this data
    fname = ['frame_',sprintf('%04d',t),'.mat'];
    if(~exist([obj.exp_info.nuc_seg_dir],'dir'))
        mkdir(obj.exp_info.nuc_seg_dir)
    end
    save([obj.exp_info.nuc_seg_dir,fname],'frame_obj')
    disp(['Writing: ',fname] )


end
    
   

end
 
%% Checking for edge effects in the contour step. 
%Basically looking to see if calculating boundary hull significantly
%increases the area. This implies contour failed at the image edge. 
function cells = check_for_edge_effects( cells, areas, img_size )

    %Loop over cells calculate boundary
    for i= 1:length(cells)
        
        cell_bnd = boundary(cells{i});
        new_cells{i} = cells{i}(cell_bnd,:);
    end

    new_areas = contour_areas( new_cells)

end
    %% Calculate the areas of all contours for filtering 
function areas = contour_areas(cs)

areas =zeros(1,length(cs));
for i = 1:length(cs)
    areas(i) = polyarea(cs{i}(:,1),cs{i}(:,2));
end

end


%% Replicate the scipy.ndimage.filters.gaussian_filter used in data_doctor
%code so that segmentation is the same here
function img = imgaussfilt_py(im, sigma )

    %Calculate a width of the gaussian based on sigma
    w = 2*floor(4*sigma + 0.5) + 1;
    
    %Now perform guassian filter
    img = imgaussfilt( im, sigma, 'FilterSize', w);
    
end



%% Remove contours within a larger contour. 
function ids = rm_duplicate_centroids(D,A,m, cells, rm_cont_dist)
    %debugging this function
    debug=0;
    %Go over distance matrix and determine if (i,j) is less than rm_dist.
    %If so, then remove it! Going over only upper triangle
    c = 0;
    rm_pair = zeros(m,m);
    idc_bad = [];
    
    %Sum over cells
    for i = 1:m
        
        %Loop over upper triangle
        for j = i+1:m
            %Pair distance
            d_pair(i,j) = D( (i-1)*(m-i/2)+j-i ); 
            if(d_pair(i,j) <= rm_cont_dist)
                rm_pair(i,j) = 1;
            end
        end
        
        %Now we can take the biggest contour if there's any overlapping
        %ones
        overlap_count = sum( rm_pair(i,:) ); %Might be multiple overlaps

        if(overlap_count > 0) 
            
            %Make sure cell(i) is not on bad list
            if( sum( idc_bad == i) == 0 )         
                %Find indices of area to query
                area_q_idx    = [i, find( rm_pair(i,:) == 1) ];
                area_comp     = A(area_q_idx);
                %Find maximum area
                [~,max_id]    = max(area_comp);
                %Now see if any of the inner cells are all within the outer
                %one
                max_cell = cells{area_q_idx(max_id)};
                if(debug)
                    figure(14)
                    hold on
                    p_blue = plot(max_cell(:,1),max_cell(:,2),'b');
                    hold off
                    label(p_blue,num2str( area_q_idx( max_id ) ))
                end
                query_cells = setdiff(area_q_idx, area_q_idx(max_id));
                for k = 1: length(query_cells)
                    if(debug)
                        hold on
                        p_red = plot(cells{ query_cells(k) }(:,1), cells{ query_cells(k) }(:,2),'r')
                        hold off
                        label(p_red,num2str( query_cells( k ) ))
                    end

                    IN = inpolygon(cells{ query_cells(k) }(:,1),cells{ query_cells(k) }(:,2), max_cell(:,1), max_cell(:,2));
                    if( sum( IN ) == length(IN) ) 
                        idc_bad = [idc_bad, query_cells(k)];
                    end
                end
                c=c+1;
                ids(c)        = area_q_idx(max_id);
            end
        else
            %No overlapping cells found, just keep i if it's not on bad
            %list
            if( sum( idc_bad == i) == 0 );
                c=c+1;
                ids(c)        = i;
            end
        end
        if(debug)
            set(gca,'Ydir','reverse')
        end
    end
    if(debug)
        ids
        idc_bad
        pause
    end
    
    %Failsafe to make sure ids aren't repeated
    ids = unique(ids);
    
end
    