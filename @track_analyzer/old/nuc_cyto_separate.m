% Calculate mean yap signal (could be arbitrary signal based on image).
% Uses the roi2poly tool. Maybe faster?
%Now computes cyto to nuclear ratio by taking a mean of a thin ring outside
%of the contour

function obj = nuc_cyto_separate(obj, params)

main_dir = obj.exp_info.sub_dir;
signal_img_dir = obj.exp_info.signal_dir;

signal_img_files =dir([signal_img_dir,'*.tif']);

track_file = [main_dir,'track_data.mat'];

debug = 0;
if(debug)
    figure(7)
    clf(7)
    %set(gcf,'color','w')    
    color_vec  = hsv(length(obj.tracks));
    color_vec  = color_vec(randperm(size(color_vec,1)),:);
end

    %For quickly querying the signal location within contour
    img = imread([signal_img_dir,signal_img_files(1).name]);
    %Generate pixel location matrices
    xgv=[1:size(img,2)];
    ygv=[1:size(img,1)];
    [X,Y]=meshgrid(xgv,ygv);


for i = 1:length(obj.img_files)
    display(['Analyzing mean signal of frame: ',num2str(i)])
    %load signal image
    img_file = [signal_img_dir,signal_img_files(i).name];
    img = imread(img_file);
    img_dbl = double(img);
    if(debug);imagesc(img);colormap('gray');end;
    
    for j = 1:length(obj.tracks)
        %if track exists at frame i
        times = obj.tracks{j}(:,1);
        if( sum( times == i) == 1 ) 
            time_idx = find(times == i);
            cell_id = obj.tracks{j}(time_idx,2);
            
            %Contour
            x     = obj.cells{i}{cell_id}(:,1);
            y     = obj.cells{i}{cell_id}(:,2);

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
            roi_bw_nuc = poly2mask(x_in,y_in,size(img,1),size(img,2));
            %Outer contour
            roi_bw_in  = poly2mask(x,y,size(img,1),size(img,2));
            roi_bw_out = poly2mask(x_out,y_out,size(img,1),size(img,2));
            %Difference contour
            roi_diff = logical(roi_bw_out - roi_bw_in);
            %Nuc/cyto ratio
            nuc = img_dbl( roi_bw_nuc );
            %sel = nuc >=500;
            mean_nuc = mean( nuc );
            mean_bg  = mean( img_dbl( roi_diff) );   
            %Save nucleus signal
            obj.tracks{j}(time_idx,5) = mean_nuc;
            %Save cytoplasm signal
            obj.tracks{j}(time_idx,6) = mean_bg; 
            
            %Debugging image calculations
%             sel_mat = zeros(size(img));
%             sel_mat( y_search(1):y_search(end),x_search(1):x_search(end) ) = IN;
%             img(logical(sel_mat)) = 500;
%             imagesc(img)

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
    end
    
    if(debug);pause;if(i < length(obj.img_files));clf(7);end;end
end

       
if(debug);hold off;end;

obj.nuc_cyto_calc = true;

end