
%Load contours, refines the contours based on knowledge of typical nucleus
%geometry. Input a directory, and this runs through. Degbug=1 for showing
%output! 

function segment_refinement_r8( exp_info, params)
% %Get data
% main_dir = '/media/jan/smbshare/Bee/Matt/cell_lines/10A_yap/zeiss_confocal/11-10-16_glass_sir-dna/CLAHE_otsu/';
% main_dir = '/Users/franklin/foobar/Bee/Matt/cell_lines/10A_yap/zeiss_confocal/11-11-16_glass_sir-dna/CLAHE_phan/';
% main_dir = '/Users/franklin/foobar/Bee/Matt/cell_lines/10A_yap/zeiss_confocal/11-11-16_glass_sir_dna/CLAHE_phan/';
% main_dir = '/media/jan/smbshare/Bee/Matt/cell_lines/10A_yap/olympus_confocal/yap_11-13-16_SirDNA_scratch/cycle_max_proj_nuc/';
% main_dir = '/media/jan/smbshare/Bee/Matt/cell_lines/10A_yap/zeiss_confocal/11-27-16_PDMS_sir_dna/nuc_CLAHE_otsu_sub/';

exp_info.nuc_seg_dir
mat_files = dir([exp_info.nuc_seg_dir,'*.tif.mat'])

%Parameters
    %lots of outputing for debugging
    %params.debug=0;
    %Plot final output
    params.final=0;
    %Negative curvature threshold
    params.curv_thresh=0.06;
    %Theta threshold
    params.theta_thresh = 40;
    %Euclidean Distance threshold
    params.dist_max_thresh = 33;
    params.dist_min_thresh = 0;
    %Max intersect dist
    params.dist_i_max= 3000
    %Distance along path threshold (number of points)
    params.path_d_min = 30;
    %If convex hull results in much greater area, use the hull 
    params.convhull_factor=1.8;
    %If a new cell has a small size, then reject!
    params.min_cell_area=200;
    %Remove contours within a min distance
    params.rm_cont_dist = 10; 
    %Removes random blobs
    params.max_circ = 2.6;
    %Minimum fraction of invagination area over A_convex - A
    params.min_invag_fraction = 0.85;
    
    %For each cell_finder data ->>>
    for f = 1:length(mat_files)

        display(['Now analyzing frame: ',num2str(f)]);

        %Load the cell finder data
        data = load([exp_info.nuc_seg_dir,mat_files(f).name]);%
        img_file = data.img;
        %Fun through the refinment operations
        run_main_fun(data,params,exp_info.nuc_seg_dir,mat_files,img_file,f)
        if(params.debug); pause; end;
        
    end

end

function run_main_fun(data,params,main_dir,mat_files,img_file,f)
    
    %Unload params
    %debug
    debug = params.debug;
    %Final 
    final = params.final;
    %Negative curvature threshold
    curv_thresh = params.curv_thresh;
    %Theta threshold
    theta_thresh=params.theta_thresh;
    %Euclidean Distance threshold
    dist_max_thresh = params.dist_max_thresh;
    dist_min_thresh = params.dist_min_thresh;
    %Max intersect dist
    dist_i_max = params.dist_i_max;
    %Distance along path threshold (number of points)
    path_d_min = params.path_d_min;
    %If convex hull results in much greater area, use the hull 
    convhull_factor = params.convhull_factor;
    %If a new cell has a small size, then reject!
    min_cell_area = params.min_cell_area;
    %Remove contours within a min distance
    rm_cont_dist = params.rm_cont_dist;     
    %max std curvature
    max_circ = params.max_circ;
    %Min area fraction of invagination
    min_invag_fraction = params.min_invag_fraction;
    
    
    img=imread(img_file);
    %Generate pixel location matrices
    xgv=[1:size(img,2)];
    ygv=[1:size(img,1)];
    [X,Y]=meshgrid(xgv,ygv);
    
    if(final)
        figure(3)
        clf(3)
    end

    if(debug)
        figure(2)
        clf(2)
        set(gcf,'color','w')
        %axis tight
    end
    
    %Quick pass to remove cells that are interior to each other
       %First calculate distance matrix based on centroids
       p=0;
       for i = 1:length(data.cells)
           p=p+1;
           centroids(p,:) = [mean( data.cells{p}(:,1)), mean( data.cells{p}(:,2) )];
           areas(p)       = polyarea( data.cells{p}(:,1), data.cells{p}(:,2) );
       end
        %centroids = cell2mat( frame_obj.centroids(:) );
        D = pdist( centroids); 
        clean_ids = rm_duplicate_centroids(D,areas,size(centroids,1),data.cells, rm_cont_dist);
        
    %Counter for accounting for newly created cells
    c=0;
    %clean_ids = 44;
    for i =  1:length(clean_ids)
        %clean_ids = 30;
        x=data.cells{clean_ids(i)}(:,1);
        y=data.cells{clean_ids(i)}(:,2);

        %Calculate centroid of the object

        %Smooth the contour
        x = smooth(x,25,'sgolay');
        y = smooth(y,25,'sgolay');

       if(debug)
            figure(2)
            hold on
            plot(x,y,'b')
            %hold on
            %Label
            L(i) = plot(mean(x),mean(y));
            label(L(i),num2str(clean_ids(i)),'fontsize',10)
            hold off
       end

        %Fix blobs with a single invagination
        %Calculate area
        A(i) = polyarea(x,y);
        Konv = flipdim(convhull(x,y),1);
        if(debug)
            display(['Blob ',num2str(i)])
        end
        
        A_conv(i) = polyarea(x(Konv),y(Konv));
            %Detect when convex hull skips over a section. Need to start at
            %2 because the convhull automatically closes shape
            area_diff = A_conv(i) - A(i);
            conv_segs = find( diff(Konv ) ~= 1);
            a_query_fraction = [];
            for j = 1: length( conv_segs ) -1 
                %Only do segs that start away from beginning of contour.
                %Sometimes there's a weird overlap at the beginning of
                %shape. 
                seg_start= Konv( conv_segs(j) );
                if(seg_start > 10)
                    seg_stop = Konv( conv_segs(j) + 1 );
                    if(debug)
                        hold on
                        plot( x(seg_start),y(seg_start),'*g')
                        plot( x(seg_stop) ,y(seg_stop) ,'*r')
                        plot( x(seg_start:seg_stop) , y(seg_start:seg_stop), 'k')
                        hold off
                    end

                    %Find area of the invagination
                    a_query = polyarea( x(seg_start:seg_stop), y(seg_start:seg_stop) );
                    %Find percent of recoverable area due to this invagination
                    a_query_fraction(j) = a_query / area_diff;
                end
            end
            
            if( max(a_query_fraction) >= min_invag_fraction )
                %This invagination makes up most of the difference. Now
                %lets use the convex hull for the whole blob. 
                x = x(Konv);
                y = y(Konv);
            end

        %Calculate local curvature (doesn't provide weather convex or concave!)
        K = LineCurvature2D([x,y]);
        
        %Remove blobs that are too curvy with a max circularity parameter
        skip_cell=0;        
        geom=polygeom(x,y);
        circularity(i) = (geom(4) .^ 2) ./ (4 * pi * geom(1));
        if( circularity(i) > max_circ )
            skip_cell=1;
            if(debug)
                hold on
                plot(x,y,'r','linewidth',3)
                hold off
            end
        end
        
        
        Vertices = [x,y];
        N = LineNormals2D(Vertices);
        
%         if(A_conv(i) >= A(i)*convhull_factor)
%             x = x(Konv);
%             y = y(Konv);
%             Vertices = [x,y];
%             N = LineNormals2D(Vertices);
%             K = LineCurvature2D([x,y]);       
%         end
        
        %Smooth curvature
        K_sm=smooth(K,3);

        %Make normal points just within or outside of polygon
        K_sign = sign(K_sm);
        K_sign(K_sign==0)=-1;
        x_norm=x + 0.5*K_sign.*N(:,1);
        y_norm=y + 0.5*K_sign.*N(:,2);

        %Figure out weather each point is convex or concave
        IN = inpolygon(x_norm,y_norm,x,y);

        %Select regions of local concavity
        sel_concav = IN;

        %Erode small features
        n_erode=2;
        sel_concav_sm = imerode(sel_concav,ones(1,n_erode)');

        %Plot smooth curvature vals on concave points
        IN=sel_concav_sm;

        %Find segments
        seg=[1; find(diff(sel_concav_sm)~=0)+1];

        %Look through segments and find center index and calc mean Normal;
        clear mean_seg_K
        clear seg_conc
        clear mean_seg_x mean_seg_y mean_N seg_median
        for j = 1:length(seg)
            seg_start=seg(j);
            if(j<length(seg))
                seg_stop = seg(j+1) - 1;
            else
                seg_stop = length(sel_concav_sm);
            end

            mean_seg_K(j) = mean( K_sm(seg_start:seg_stop) );
            %Weight x,y,N by curvature? 
            k_vals = K_sm(seg_start:seg_stop);
            weights = k_vals ./ sum(k_vals);
            seg_median(j) = round( mean([seg_start,seg_stop]));
            mean_seg_x(j) = sum(x(seg_start:seg_stop).*weights);
            mean_seg_y(j) = sum( y(seg_start:seg_stop).*weights);
            mean_N(j,:)   = sum( N(seg_start:seg_stop,:).*[weights,weights],1);
            %Determine if seg is concave
            seg_conc(j)   = sum(sel_concav_sm(seg_start:seg_stop)) == length(k_vals);

            if(abs(mean_seg_K(j)) > curv_thresh & seg_conc(j)==1)
                if(debug)
                    hold on
                    plot([mean_seg_x(j) mean_seg_x(j)+mean_seg_K(j).*10*mean_N(j,1)']',[mean_seg_y(j) mean_seg_y(j)+mean_seg_K(j).*10*mean_N(j,2)']','k');
                    hold off
                end
            end
        end

        %Use segments above threshold
        seg_sel_idx = find(mean_seg_K > curv_thresh & seg_conc == 1);
        clear theta theta_sel
        clear path_dist path_dist_sel
        clear sel_mat
        clear dist dist_sel
        clear vec_oppo
        clear dist_i_sel dist_i

        if(sum(seg_sel_idx)>0 && skip_cell==0)
            %Run through selected segments and see if they are opposing within
            %threshold angle? 
            %mean_N(seg_sel_idx,:) = abs(mean_N(seg_sel_idx,:));
            for j = 1:length(seg_sel_idx)
                if(debug);hold on; h(j) = plot([mean_seg_x(seg_sel_idx(j)) mean_seg_x(seg_sel_idx(j))+mean_seg_K(seg_sel_idx(j)).*50*mean_N(seg_sel_idx(j),1)']',[mean_seg_y(seg_sel_idx(j)) mean_seg_y(seg_sel_idx(j))+mean_seg_K(seg_sel_idx(j)).*50*mean_N(seg_sel_idx(j),2)']','k');
                label(h(j),num2str(j));hold off;end %<<< for diagnostic purposes
                vec_a = mean_N(seg_sel_idx(j),:);
                for k = 1:length(seg_sel_idx)
                    %Calculate intersection point between vectors (basically
                    %accounts for theta and direction
                    x1 = mean_seg_x(seg_sel_idx(j));
                    x2 = x1 + mean_seg_K(seg_sel_idx(j))*mean_N(seg_sel_idx(j),1);
                    y1 = mean_seg_y(seg_sel_idx(j));
                    y2 = y1 + mean_seg_K(seg_sel_idx(j))*mean_N(seg_sel_idx(j),2);

                    x3 = mean_seg_x(seg_sel_idx(k));
                    x4 = x3 + mean_seg_K(seg_sel_idx(k))*mean_N(seg_sel_idx(k),1);
                    y3 = mean_seg_y(seg_sel_idx(k));
                    y4 = y3 + mean_seg_K(seg_sel_idx(k))*mean_N(seg_sel_idx(k),2);                

                    %Find intersection point
                    mA = (y1-y2)/(x1-x2);
                    mB = (y3-y4)/(x3-x4);
                    ua = ( (x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)  ) ./( (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1) );
                    x_i = x1 + ua*(x2 - x1);
                    y_i = y1 + ua*(y2 - y1);

                    [x1,y1; x2, y2; x3, y3; x4,y4];

                    %Calculate mean_distance
                    dist_i(j,k) = mean( [(x1-x_i)^2+(y1-y_i)^2, (x3-x_i)^2+(y3-y_i)^2]);

                    %Check by plotting
                    %plot([x1,x2],[mA*(x1-x1)+y1,mA*(x2-x1)+y1])
                    %plot([x3,x4],[y3,y4])
                    %Check by plotting
                    %plot(x_i,y_i,'*','markersize',10)

                    vec_b = mean_N(seg_sel_idx(k),:);
                    %Calculate angle between vector j and all others 
                    a_norm         = vec_a ./ mag(vec_a);
                    b_norm         = vec_b ./ mag(vec_b);
                    theta(j,k)     = acos(sum( a_norm.*b_norm ));
                    dist(j,k)      = sqrt( (mean_seg_x(seg_sel_idx(j)) - mean_seg_x(seg_sel_idx(k)))^2 + (mean_seg_y(seg_sel_idx(j)) - mean_seg_y(seg_sel_idx(k)))^2);

                    %Calculate distance along path
                    path_dist(j,k) = abs( seg(seg_sel_idx(j)) - seg(seg_sel_idx(k)));
                end
            end

            %Create selection matrices   
            theta  = theta .* 180 ./ pi;
            dist_i_sel = dist_i <= dist_i_max;
            path_dist_sel = path_dist >= path_d_min;
            theta_sel = theta <= theta_thresh | theta >= 180 - theta_thresh;
            dist_sel  = dist < dist_max_thresh & dist > dist_min_thresh;

            %If more than two points, then try dist_i_sel
            if(size(theta,1)>2)
                sel_mat = theta_sel & dist_sel & path_dist_sel & dist_i_sel;
            else
                sel_mat = theta_sel & dist_sel & path_dist_sel ;
            end

            %output for debugging
            if(debug)
                theta
                dist_i
                dist
                dist_sel
                theta_sel
                dist_i_sel
                path_dist_sel
                sel_mat
            end        

            %Determine pairs, if there are pairs that meet criteria, then split
            %cells
            if(sum(sel_mat(:))>0 & sum(sel_mat(:))<=8)
                pairs  = find_pairs(sel_mat);

                %Pair validator - check if index is used more than once! 
                pairs = check_pairs(pairs, dist, theta,path_dist);

                %For each pair, pinch cell off at 'seg_median' 
                n_new_cells=size(pairs,1)+1;
                %First make cells out of the continuous stretches (outer cells - always 2)
                if(n_new_cells>2)
                    for j = 1:n_new_cells - 1 %2 %size(pairs,1)
                        median_val_a = seg_median(seg_sel_idx(pairs(j,1)));
                        median_val_b = seg_median(seg_sel_idx(pairs(j,2)));
                        %Check if these are right
                        %plot(x(median_val_a),y(median_val_a),'*r','markersize',10);
                        %plot(x(median_val_b),y(median_val_b),'*r','markersize',10);
    
                        %Start with pair #1. Find the seg start point.
                        seg_start=min( seg_median(seg_sel_idx(pairs(j,:))) );
                        %Guess the stop point. 
                        seg_stop =max( seg_median(seg_sel_idx(pairs(j,:))) );
                        %Check to see if the other points are between the
                        %seg_start and seg_stop
                            seg_indices = [seg_start, seg_stop];
                            p=0;
                            for k = 1:size(pairs,1)
                                if(k~=j)
                                    p=p+1;
                                    seg_indices = [seg_indices, seg_median(seg_sel_idx(pairs(k,:)))];
                                end
                            end
                        seg_indices = sort(seg_indices);
                        tmp_sel = seg_indices > seg_start & seg_indices < seg_stop;
                        cells_between = sum(tmp_sel) / 2;
                        %If there are cells between, then make the
                        %appropriate segments
                        seg_ids = [];
                        if(cells_between > 0)
                            seg_ids(1) = seg_start;
                            tmp_id     = find(seg_indices == seg_start);
                            seg_ids(2) = seg_indices(tmp_id + 1);
                            tmp_id     = find(seg_indices == seg_stop);
                            seg_ids(3) = seg_indices(tmp_id - 1);
                            seg_ids(4) = seg_stop;
                            
                            %Build cell
                            data_tmp = [x(seg_ids(1):seg_ids(2)),y(seg_ids(1):seg_ids(2))];
                            data_tmp = [data_tmp; x(seg_ids(3):seg_ids(4)),y(seg_ids(3):seg_ids(4))];
                            
                        else
                            data_tmp = [x(seg_start:seg_stop),y(seg_start:seg_stop)];
                        end
                        
                        %Check size
                        a_tmp = polyarea( data_tmp(:,1), data_tmp(:,2) );
                        if(a_tmp > min_cell_area)
                        %Copy contour data
                            c=c+1;
                            data_refined.cells{c} = data_tmp;
                            data_refined.area{c} = a_tmp;
                            %Plot contour
                            data_refined.cells{c}(:,:);
                            if(debug)
                                hold on; 
                                plot(data_refined.cells{c}(:,1),data_refined.cells{c}(:,2),'--r','linewidth',4);
                                hold off
                            end

                        end
                    end
                    %Now connect the non-continuous cell (only one cell!) 
                        %Start with end of pair one and beginning of pair two
                        %seg_indices = sort(seg_median(seg_sel_idx(pairs(:))));
                        
                        seg_ids = [];
                        seg_ids(1) = 1;
                        seg_ids(2) = seg_median(seg_sel_idx(1));
                        seg_ids(3) = seg_median(seg_sel_idx(end));
                        seg_ids(4) = length(x);
                        
                        data_tmp = [x(seg_ids(1):seg_ids(2)),y(seg_ids(1):seg_ids(2))];
                        data_tmp = [data_tmp; x(seg_ids(3):seg_ids(4)),y(seg_ids(3):seg_ids(4))];
                        
                            %Check size
                            a_tmp = polyarea(data_tmp(:,1),data_tmp(:,2));
                            if(a_tmp > min_cell_area)                        
                                c=c+1;
                                data_refined.cells{c}=data_tmp;
                                data_refined.area{c}=a_tmp;
                                if(debug)
                                    hold on
                                    plot(data_refined.cells{c}(:,1),data_refined.cells{c}(:,2),'--b','linewidth',5);
                                    hold off
                                end
                            end
                else
                    %Case with only 2 new cells
                    %Start with end of pair one and beginning of pair two
                    seg_indices = sort(seg_median(seg_sel_idx(pairs(:))));
                    %First cell (continuous section)
                        %Check area
                        a_tmp = polyarea(x(seg_indices(1):seg_indices(2)),y(seg_indices(1):seg_indices(2)));
                        if(a_tmp > min_cell_area)
                            c=c+1;
                            data_refined.cells{c} = [x(seg_indices(1):seg_indices(2)),y(seg_indices(1):seg_indices(2))];
                            data_refined.area{c} = a_tmp;
                            if(debug)
                                hold on
                                plot(data_refined.cells{c}(:,1),data_refined.cells{c}(:,2),'--r','linewidth',5);
                                hold off
                            end
                        end

                    %Second cell (broken section)
                        data_tmp = [[x(seg_indices(2)+1:end),y(seg_indices(2)+1:end)]];
                        data_tmp = [data_tmp; [x(1:seg_indices(1)-1),y(1:seg_indices(1)-1)]];
                        %Check area
                        a_tmp = polyarea(data_tmp(:,1),data_tmp(:,2));
                        if(a_tmp > min_cell_area)
                            c=c+1;
                            data_refined.cells{c}=data_tmp;
                            data_refined.area{c} = a_tmp;
                            if(debug)
                                hold on
                                plot(data_refined.cells{c}(:,1),data_refined.cells{c}(:,2),'--b','linewidth',5);
                                hold off
                            end
                        end
                 end
            else
                c=c+1;
                data_refined.cells{c}=[x,y];
                data_refined.area{c} = A(i);
            end
        elseif(skip_cell==0)
            c=c+1;
            data_refined.cells{c}=[x,y];
            data_refined.area{c} = A(i);
        end
    end

    %Plot centered segment norms
    %plot([Vertices(sel_concav_sm,1) Vertices(sel_concav_sm,1)+K_sm(sel_concav_sm).*N(sel_concav_sm,1)]',[Vertices(sel_concav_sm,2) Vertices(sel_concav_sm,2)+K_sm(sel_concav_sm).*N(sel_concav_sm,2)]','r');    
    hold off
    axis equal   
    set(gca,'ydir','reverse')

%     %Last step is to remove contours within other contours!
%         %First calculate distance matrix based on centroids
%         centroids = cell2mat( frame_obj.centroids(:) );
%         areas     = cell2mat( frame_obj.area );          
%         D = pdist( centroids); 
%         clean_ids = rm_duplicate_centroids(D,areas,size(centroids,1),rm_cont_dist);

    %Remove contours again in case the cell splitting algorithm created new
    %cells with interior blobs
    p=0;
    centroids=[];
    areas   =[];
    for i = 1:length(data_refined.cells)
       p=p+1;
       centroids(p,:) = [mean( data_refined.cells{p}(:,1)), mean( data_refined.cells{p}(:,2) )];
       areas(p)       = polyarea( data_refined.cells{p}(:,1), data_refined.cells{p}(:,2) );
    end
    %centroids = cell2mat( frame_obj.centroids(:) );
    D = pdist( centroids); 
    clean_ids = rm_duplicate_centroids(D,areas,size(centroids,1),data_refined.cells, rm_cont_dist);
  
    if(final); figure(3); end;

    k=0;
    for i = 1:length(clean_ids)
        k=k+1;
        %Calculate centroid
        frame_obj.centroids{k} = [mean(data_refined.cells{clean_ids(i)}(:,1)),mean(data_refined.cells{clean_ids(i)}(:,2))];
        %Add area to frame_obj
        frame_obj.area  = data_refined.area;

        if(final)
            hold on
            h(i) = plot(data_refined.cells{clean_ids(i)}(:,1), data_refined.cells{clean_ids(i)}(:,2),'r');
            %label(h(i),num2str(i));
            hold off
        end
    end

    %Append the new cell data
    frame_obj.refined_cells = data_refined.cells(clean_ids);        
    
    if(final); set(gca,'Ydir','reverse');end;
    frame_obj.img_file      = img_file;
    save([main_dir,mat_files(f).name],'frame_obj','-append');
end

%Remove contours within a larger contours. 
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
end

%Find next closest indices
function N = next_closest(q,indices)
    indices;
    sel = find(indices > q);
    diff_array = indices(sel) - q;
    [~,n] = min(diff_array);
    N = indices(sel(n));
end

%Magnitude
function M = mag(v)
    M = sqrt(sum(v.^2));
end

%Check pairs of points for reused points. Then take the pair with the
%minimum distance
function new_pairs = check_pairs(pairs, dist, theta,path_dist)
    %Output
    new_pairs = [];
    %Previously checked idx
    skip=[];
    idx = pairs(:);
    tmp=0;
    for i = 1:length(idx)
        %tmp = tmp + sum(sum(pairs == idx_2_chk(i)));
        %display(['Checking point ',num2str(idx(i))])
        [r,c] = find(pairs == idx(i));
        if(length(r) > 1 )

            %skip=[skip;pairs(r(1),c(1))];
            %Counter for repeats
            tmp=tmp+1;  
            %Computer distance between pairs
            dist_metric = [];
            row_vec     = [];
            k=0;
            for j= 1:length(r)
                %Check if any points in the pairs match skip
                pairs(r(j),:);
                chk = sum(skip==pairs(r(j),1)) + sum(skip==pairs(r(j),2));
                if(chk == 0);
                    k=k+1;
                    query_idx = pairs(r(j),:);
                    dist_metric(k) = dist( query_idx(1),query_idx(2) );
                    row_vec(k)     = r(j);
                end
            end
            %Select best pair based on closest distance between
            %repeated point
            if( ~isempty(dist_metric))
                [~,sel] = min(dist_metric);            
                new_pairs = [new_pairs; pairs(row_vec(sel),:)];
                %Now also add the point paired to idx(i) to the skip list
                skip = [ skip, pairs(row_vec(sel),:) ] ;
            end
            
        end

    end
    
    if(tmp == 0)
        new_pairs = pairs;
    end
    
end


%Determine pairs of points for pinching cell off
function pairs = find_pairs(sel_mat)
    [x,y] = find(sel_mat==1);
    pairs = zeros(length(x)/2,2);
    k=0;
    for i = 1:length(x)
        pair_q = [x(i),y(i)];
        for j = 1:size(pairs,1)
            chk1 = pair_q(1) ~= pairs(j,1) | pair_q(2) ~= pairs(j,2);
            chk2 = pair_q(1) ~= pairs(j,2) | pair_q(2) ~= pairs(j,1); 
            if( chk1 & chk2 )
                chk(j) = 0;
            else
                chk(j) = 1;
            end
        end
        if( sum(chk) == 0)
            %Keep pair
            k=k+1;
            pairs(k,:)=pair_q;
        end
    end
end     
