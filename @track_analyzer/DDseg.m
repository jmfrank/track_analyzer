%Nucleus segmentation automated for track_analyzer

function DDseg( obj, extra )


%Get exp info
exp_info = obj.exp_info;

%Get list of DD seg files
seg_files = dir([exp_info.nuc_seg_dir,'*.tif.mat']);

    
%Now loop over these files and segment
for i = 1:length(seg_files)
    
    %Clear variables?
    clear frame_obj
    
    disp( seg_files(i).name )
    
    %Loading datafile
    data = load([exp_info.nuc_seg_dir,seg_files(i).name]);
        
    %load img 
    img = imread(data.img);
    
    
    %If extra is specified, do some extra refinement steps
    if( exist('extra','var') )

        %Quick pass to remove cells that are interior to each other
        %First calculate distance matrix based on centroids
        clear centroids areas
        for p = 1:length(data.cells)
           centroids(p,:) = [mean( data.cells{p}(:,1)), mean( data.cells{p}(:,2) )];
           areas(p)       = polyarea( data.cells{p}(:,1), data.cells{p}(:,2) );
        end
        D = pdist( centroids); 
        clean_ids = rm_duplicate_centroids(D,areas,size(centroids,1),data.cells, extra.rm_cont_dist);



        %Create empty image
        bw = zeros( size(img,1),size(img,2));

        %Counter
        c = 0;

       %Loop over potential cells
       for j = clean_ids
          
           %temporary cell
           tmp_cell = data.cells{j};
           
%            %Dilate tmp_cell if expand_ratio specified
%            if( isfield(extra,'expand_ratio') )
%                %Extract data
%                x = data.cells{j}(:,1);
%                y = data.cells{j}(:,2);
%                
%                %Transform
%                A = extra.expand_ratio/100 + 1;
%                
%                x = A*x+(1-A)*mean(x);
%                y = A*y+(1-A)*mean(y);
%                %Replace
%                tmp_cell= [x,y];
%            end
           
           %Get mean intensity of this cell from raw image
           tmp_bw = poly2mask( data.cells{j}(:,1),data.cells{j}(:,2),size(img,1),size(img,2));
           tmp_int = mean( img( tmp_bw ) ); %<< mean intensity of this cell
           
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
    [a,b,c] = fileparts([exp_info.nuc_seg_dir,seg_files(i).name]);
    frame_file = [a,'/',b,'.mat'];
    save(frame_file,'frame_obj','-append')
end

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
                ids(c) = i;
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
