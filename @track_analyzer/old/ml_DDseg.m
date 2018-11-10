%Simple version of data_doctor cell segmentation. 

function obj = ml_DDseg(obj, params, extra)


debug = 0;

Z = obj.exp_info.z_frames;
T = obj.exp_info.t_frames;

%Get filters
    %Preguassian
    pre_gauss  = fspecial('gaussian',params.pre_gauss);
    %Dog fiter
    dog_filter_big = fspecial('gaussian', params.dog(1));
    dog_filter_small = fspecial('gaussian', params.dog(2));

    

%LSM time series goes Z first, then time. For each time point, load the
%Z-stack at time t, then max project and segment. 

for t = 1:T
    
    %Clear out frame_obj
    frame_obj = struct('refined_cells',[],'centroids',[],'BW',[]);

    %Calculate indices. Skip every other b/c thumbnails? 
    start = Z.*(t-1)*2 + 1;
    indices  = start:2:start + 2*Z-1;
    %Grab this Z-stack
    IMG = tiffread( obj.exp_info.lsm_file, indices );
    
    [y,x] = size( IMG(1).data );
    %Now project
    img = zeros(y,x,Z);
    for z = 1:Z
        img(:,:,z) = IMG(z).data;
    end
    mip = max(img,[],3);
 
    %Prefilter
    mip_pre = imgaussfilt_py( mip, params.pre_gauss);
    
    %Filter mip image
    mip_dog = imgaussfilt_py( mip_pre,params.dog(1)) - imgaussfilt_py(mip_pre, params.dog(2));
    
    %Get contours
    C = contourc(mip_dog, [params.dog_thresh,params.dog_thresh]);
    cells = C2xyz(C);
    %Get areas 
    areas = contour_areas( cells );
    
    %Filter areas 
    sel = areas >= params.min_area & areas <= params.max_area;
    cells = cells(sel);
    areas = areas(sel);
    
    %If extra is specified, do some extra refinement steps
    if( exist('extra','var') )

        %Quick pass to remove cells that are interior to each other
        %First calculate distance matrix based on centroids
        centroids = cell2mat(cellfun(@(x) mean(x),cells,'UniformOutput',false)');
        D = pdist( centroids); 
        clean_ids = rm_duplicate_centroids(D,areas,size(centroids,1),cells, extra.rm_cont_dist);

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
           
           %Use threshold to ensure this is a fluorescence positive cell...
           if tmp_int >= extra.int_thresh 
            c=c+1;
            frame_obj.refined_cells{c} = tmp_cell;
            frame_obj.centroids{c} = mean( tmp_cell );
            bw = bw + tmp_bw;
           end
           
       end
       
        
    else %Just use all data from DD
        
        frame_obj.refined_cells = data.cells;
        bw = zeros( size(img,1),size(img,2));
        for j = 1:length(data.cells)
            %Just take all data from datadoc
            frame_obj.centroids{j} = mean( data.cells{j});
            bw = bw + poly2mask(data.cells{j}(:,1),data.cells{j}(:,2),size(img,1),size(img,2));
        end
    end
    %Add the BW image to frame_obj
    frame_obj.BW = bw;

    %Now save this data
    fname = ['frame_',sprintf('%04d',t),'.mat'];
    if(~exist([obj.exp_info.nuc_seg_dir],'dir'))
        mkdir(obj.exp_info.nuc_seg_dir)
    end
    save([obj.exp_info.nuc_seg_dir,fname],'frame_obj')
    disp(['Writing: ',fname] )

    %Display in debug
    if(debug)
        subplot(1,3,1)
        imagesc(mip_dog)
        colormap gray
        subplot(1,3,2)
        contour(mip_dog,[params.dog_thresh,params.dog_thresh]);
        set(gca,'ydir','reverse')
        subplot(1,3,3)
        imagesc( bw );
        pause
    end    
    
end

end


%Calculate the areas of all contours for filtering 
function areas = contour_areas(cs)

areas =zeros(1,length(cs));
for i = 1:length(cs);
    areas(i) = polyarea(cs{i}(:,1),cs{i}(:,2));
end

end


%Replicate the scipy.ndimage.filters.gaussian_filter used in data_doctor
%code so that segmentation is the same here
function img = imgaussfilt_py(im, sigma )

    %Calculate a width of the gaussian based on sigma
    w = 2*floor(4*sigma + 0.5) + 1;
    
    %Now perform guassian filter
    img = imgaussfilt( im, sigma, 'FilterSize', w);
    
end



%Remove contours within a larger contour. 
function ids = rm_duplicate_centroids(D,A,m, cells, rm_cont_dist)
    %debugging this function
    debug=0;
    %Go over distance matrix and determine if (i,j) is less than rm_dist.
    %If so, then remove it! Going over only upper triangle
    c = 0;
    rm_pair = zeros(m,m);
    idc_bad = [];
    for i = 1:m
        for j = i+1:m
            %Pair distance
            d_pair(i,j) = D( (i-1)*(m-i/2)+j-i ); 
            if(d_pair(i,j) <= rm_cont_dist)
                rm_pair(i,j) = 1;
            end
        end
        
        %Now we can take the biggest contour if there's any overlapping
        %ones
        overlap_count = sum( rm_pair(i,:) );

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