% Calculate mean yap signal (could be arbitrary signal based on image).
% Uses the roi2poly tool. Maybe faster?
%Now computes cyto to nuclear ratio by taking a mean of a thin ring outside
%of the contour. 
%This version calculates nuc/cyto ratio for yap and tead channels. 

function obj = nuc_cyto_tead_yap(obj, params)

main_dir = obj.exp_info.sub_dir;

signal_img_1_files =dir([obj.exp_info.signal_dir,'*.tif']);
signal_img_2_files =dir([obj.exp_info.nuc_seg_dir,'*.tif']);

track_file = [main_dir,'track_data.mat'];

debug = 0;
if(debug)
    figure(7)
    clf(7)
    %set(gcf,'color','w')    
    color_vec  = hsv(length(obj.tracks));
    color_vec  = color_vec(randperm(size(color_vec,1)),:);
end


for i = 1:length(signal_img_1_files)
    display(['Analyzing mean signal of frame: ',num2str(i)])
    %load signal image 1
    img_file_1 = [obj.exp_info.signal_dir,signal_img_1_files(i).name];
    img_1 = imread(img_file_1);
    img_dbl_1 = double(img_1);
    %load signal image 2
    img_file_2 = [obj.exp_info.nuc_seg_dir,signal_img_2_files(i).name];
    img_2 = imread(img_file_2);
    img_dbl_2 = double(img_2);
    

    
    if(debug);imagesc(img_1);colormap('gray');end;
    
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
            roi_bw_nuc = poly2mask(x_in,y_in,size(img_1,1),size(img_1,2));
            %Outer contour
            roi_bw_in  = poly2mask(x,y,size(img_1,1),size(img_1,2));
            roi_bw_out = poly2mask(x_out,y_out,size(img_1,1),size(img_1,2));
            %Difference contour
            roi_diff = logical(roi_bw_out - roi_bw_in);
            %First measure YAP channel ('signal')
                mean_nuc = mean( img_dbl_1( roi_bw_in ) );
                mean_bg  = mean( img_dbl_1( roi_diff) );   
                %Write nucleus and subtract background signal
                obj.tracks{j}(time_idx,5) = (mean_nuc - params.signal_channel_bg);
                obj.tracks{j}(time_idx,6) = (mean_bg  - params.signal_channel_bg); 
            %Now measure the TEAD channel ('nuclear')
                mean_nuc = mean( img_dbl_2( roi_bw_in ) );
                mean_bg  = mean( img_dbl_2( roi_diff) );   
                %Write nucleus and subtract background signal
                obj.tracks{j}(time_idx,7) = (mean_nuc - params.nuclear_channel_bg);
                obj.tracks{j}(time_idx,8) = (mean_bg  - params.nuclear_channel_bg);

            %Plot during debugging
            if(debug)
                figure(12)
                subplot(1,3,1)
                imshow(roi_bw_nuc);
                subplot(1,3,2)
                imshow(roi_diff);
                subplot(1,3,3)
                imagesc(img_1);
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