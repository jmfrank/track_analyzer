% 3D Segmentation of cell nuclei. Based on: 'Object segmentation and ground
% truth in 3D embryonic imaging'. 
% 4-01-21: Now just a single script for 2D and 3D. 
% 3-09-21: Now object must contain steps/params for segmentation. 
% 7-23-18: changed the filtering step. Using diffusion in 3D code. Also,
% calculating gaussian gradient in 3D using separable filters. Now we don't
% have to loop over z. Hopefully this is a bit faster! ALSO, changing the
% thresholding step to use an absolute value rather than fraction of total
% luminance. 
% Also, keeping reader open during loop over time. Might be taking a long
% time to open and close image reader. 

function obj = segment(obj, seg_info, debug)


%% Pre-processing using bio-formats. 

%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = get_memo_reader(obj.exp_info.img_file);
series = 1; % Need to further define if multiple series in bf format. 

%Get the image size of this series. 
T = reader.getSizeT;
%Experiment info to pass along. 
exp_info = obj.exp_info;

% Figure out what channel and 'type' of segmentation. 
if nargin < 2
    
    debug=0;
    seg_info.type='cells';

    % Determine which channel to perform cell segmentation on. 
    chs = fieldnames(exp_info.steps);
    seg_info.channel = [];
    for i = 1:length(chs)
        name=chs{i};
        if isfield(exp_info.steps.(name),'cells')
            seg_info.channel  =str2num(name(9:end));
            break
        end
    end

elseif nargin < 3
    debug = 0;
end
       
% Write out which channel is being analyzed. 
disp(['Using channel ', num2str(seg_info.channel)])
if isempty(seg_info.channel)
    errordlg('No cell segmentation channel found')
end

%Now run loops. Params is passed around in case things change from user. 
disp(['Segmenting ',num2str(T),' frames.'])
for t = 1:T
   inner_function( obj, seg_info, t, reader, debug);
end

%Close reader at very end. 
reader.close();

%Always clear out flags in after new segmentation. 
obj = obj.clear_flags;
obj.save; %Auto-save. 

end

%% Inner function to run 
function inner_function( obj, seg_info, t, reader, debug)
    
%Display 
disp(['Started frame: ',num2str(t)]);      

%%%%%%%%%%%%%%% START PROCESSING %%%%%%%%%%%%%%%%
    image_bits     = reader.getBitsPerPixel;
    CHANNEL_name   = ['channel_', pad(num2str(seg_info.channel), 2, 'left', '0')];
    seg_type       = seg_info.type;
    %Get the bio-formats image index corresponding to this z-stack:
    I = get_stack(reader, t, seg_info.channel);
    original_size=size(I);
    %Reduce i for debugging
    if(debug)
        I = I(1:512,1:512,:);
    end

    disp( ['Segment dimension: ',num2str(obj.exp_info.params.(CHANNEL_name).(seg_type).seg_dims)] )
    
    % Check seg_dims. Max project if 2D. 
    if obj.exp_info.params.(CHANNEL_name).(seg_type).seg_dims == 2
        I = max(I,[],3);
    end
    
    %Check if we reduce image scale. 
    img_scale = obj.exp_info.params.(CHANNEL_name).(seg_type).img_scale;
    if img_scale > 1
        I = imresize3(I,1/img_scale);
        display(["Rescaling image by: ", string(img_scale)]);
        display(["Rescaling image to: ", string(size(I))]);
        
    end


    % Get mask, and initialize segmenter object. Input to segmenter must
    % have the final scaling. However, we have a function to reverse
    % scaling for output. 
    if any( strcmp('add_mask',{obj.exp_info.steps.(CHANNEL_name).(seg_type)(:).type}) )
        % Get Mask image. 
        F = obj.get_frame_files;
        load(F{1})
        % get mask channel from steps. 
        idx = find(strcmp('add_mask',{obj.exp_info.steps.(CHANNEL_name).(seg_type)(:).type}),1); 
        mask_channel = obj.exp_info.steps.(CHANNEL_name).(seg_type)(idx).value;    
        msk_channel_str = ['seg_channel_',pad(num2str(mask_channel),2,'left','0')];
        % No longer descaling masks to original size. that confused
        % everything. 
        %if img_scale > 1

            % Make a new BW from original mask data. For now, masks are
            % always cells. 
            %og_mask = cat(1,frame_obj.(msk_channel_str).cells.stats.PixelIdxList);
            %og_bw = obj.rebuild_BW(og_mask);
            %newBW=imresize3(og_bw, 1/img_scale);
            %mask_stats=regionprops(newBW,'PixelIdxList');
      
        %else
        %end


        n_cells = length(frame_obj.(msk_channel_str).cells.stats); 
        mask_stats(n_cells,1) = struct('PixelIdxList',[]);
        [mask_stats.PixelIdxList] = deal(frame_obj.(msk_channel_str).cells.stats.PixelIdxList);
    
        % send to segmenter. convert mask stats to a cell array for use
        % with segmenter. 
        mask = list_2_cell_array( mask_stats, 'PixelIdxList');
        
        S = segmenter(I, image_bits, obj.exp_info.params.(CHANNEL_name).(seg_type), t, seg_type, mask);
    else
        S = segmenter(I, image_bits, obj.exp_info.params.(CHANNEL_name).(seg_type), t, seg_type);
    end
    
    % Loop over steps. params from obj.obj.exp_info. 
    seg_steps = obj.exp_info.steps.(CHANNEL_name).(seg_type);
    
    % ignore 'add mask' step. Done above ^^^
    ignore = strcmp('add_mask',{obj.exp_info.steps.(CHANNEL_name).(seg_type)(:).type});
    seg_steps = seg_steps(~ignore);
    
    % loop over steps. 
    for i = 1:length(seg_steps)

        disp(seg_steps(i).type);
        S.(seg_steps(i).type);
    end

    % Compensate if original image was re-scaled for speed. 
    %if img_scale > 1
    %    S.reverse_scaling(original_size);
    %    display('reverse_scaled')
    %end
    

    % Assign final objects to cell masks. First we must add back the
    % original cell mask prior to scaling. 
    if strcmp(seg_type, 'spots')
        S.assign_objects_2_cells();
    end


    %% Output. 
    % Create frame_obj and save. 
    seg_channel_name=strcat('seg_',CHANNEL_name);
    frame_obj = struct(seg_channel_name,[]);
    
    % change post-processing depending on seg_type. 
    switch seg_type
        
        case 'cells'

            frame_obj.(seg_channel_name).cells.touches_border = false([1,length(S.stats)]);

          
            borders = border_frame( size(S.BW) );
            borders_idx = find(borders>0);
            
            % saving stats object instead of re-writing pixelidxlist? 
            frame_obj.(seg_channel_name).cells.stats = S.stats;

            for i = 1:length(S.stats)       
                
                %frame_obj.(seg_channel_name).cells.PixelIdxList{i} = S.stats(i).PixelIdxList;
                %do any pixels of this cell overlap with border pixels. 
                overlap = intersect(borders_idx,frame_obj.(seg_channel_name).cells.stats(i).PixelIdxList);
                if ~isempty(overlap)
                    frame_obj.(seg_channel_name).cells.touches_border(i) = 1;
                end
            end

            display('identified border touching objects')
            %Add final binarized image to frame_obj for save keeping
            %frame_obj.(seg_channel_name).cells.BW = S.BW;

            % %Trace boundaries.
            % cBW = max(S.BW, [], 3);
            % C = bwboundaries(cBW);
            % if ~isempty(C)
            %     C = cellfun(@(x) [smooth(x(:,2)),smooth(x(:,1))],C,'uniformoutput',0);
            %     c_ctr = cellfun(@(x) [mean(x(:,1)),mean(x(:,2))],C,'uniformoutput',0);
            %     %Match centroids. 
            %     fobj_ctrs = cat(1,frame_obj.(seg_channel_name).cells.centroids{:});
            %     D = pdist2(cat(1,c_ctr{:}),fobj_ctrs(:,1:2));
            %     [~,idx] = min(D);
            %     frame_obj.(seg_channel_name).cells.contours=C( idx );
            % end
            

        case 'foci'
        
            %Add final binarized image to frame_obj for save keeping
            %frame_obj.(seg_channel_name).foci.BW = S.BW;
            %Add sub-nuclear segmented nucleoli. 
            frame_obj.(seg_channel_name).foci.stats =S.stats;
            if ~isempty(S.spot_fits)
                frame_obj.(seg_channel_name).foci.fits =S.spot_fits;
            end
        case 'spots'

            frame_obj.(seg_channel_name).spots.stats=S.stats;
            if ~isempty(S.spot_fits)
                frame_obj.(seg_channel_name).spots.fits =S.spot_fits;
            end
    end
       
    %Save frame_obj
    fname = ['frame_',sprintf('%04d',t),'.mat'];
    
    if(~exist(obj.exp_info.nuc_seg_dir,'dir'))
        mkdir(obj.exp_info.nuc_seg_dir)
    end
        
    parsave([obj.exp_info.nuc_seg_dir,fname],frame_obj,seg_channel_name,seg_type)
    if exist('disp_str','var'); clearString(disp_str); end
    disp_str = ['Finished frame:     ',num2str(t)];
    disp(disp_str)
        
end

function plot_text(stats)

for i = 1:length(stats)
   text(stats(i).Centroid(1),stats(i).Centroid(2),num2str(i),'color','w') 
end

end

%Parallel saving technique
function parsave(fname, frame_obj, channel_str,seg_type)

if exist( fname, 'file')
    try
        D = load(fname,'frame_obj');
    catch
        D = struc();
        warning('encountered empty frame file')
    end
    
    %Append the newly segmented data to existing frame. 
    D.frame_obj.(channel_str).(seg_type) = frame_obj.(channel_str).(seg_type);
    
    frame_obj = D.frame_obj;

    save(fname, 'frame_obj','-v7.3');
else
    save(fname, 'frame_obj','-v7.3');

end

end



%% Plotting text over image
function    plot_liar_text( stats, liar_list )

delete(findobj(gca,'Type','text'))
hold on
for i = 1:length(liar_list)
    ctr = stats( liar_list(i) ).Centroid;
    text(ctr(1),ctr(2),num2str( liar_list(i)),'color','g')
end
hold off

end