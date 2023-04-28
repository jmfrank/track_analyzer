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
reader = bfGetReader(obj.exp_info.img_file);
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
    if isempty(seg_info.channel)
        errordlg('No cell segmentation channel found')
    end
elseif nargin < 3
    debug = 0;
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
    
    %Reduce i for debugging
    if(debug)
        I = I(1:512,1:512,:);
    end

    disp( ['Segment dimension: ',num2str(obj.exp_info.params.(CHANNEL_name).(seg_type).seg_dims)] )
    
    % Check seg_dims. Max project if 2D. 
    if obj.exp_info.params.(CHANNEL_name).(seg_type).seg_dims == 2
        I = max(I,[],3);
    end
    
    % Get mask, and initialize segmenter object. 
    if obj.exp_info.steps.(CHANNEL_name).mask_channel>0
        % Get Mask image. 
        F = obj.get_frame_files;
        load(F{1})
        msk_channel_str = ['seg_channel_',pad(num2str(obj.exp_info.steps.(CHANNEL_name).mask_channel),2,'left','0')];
        mask = frame_obj.(msk_channel_str).PixelIdxList;
        mask_centroids = frame_obj.(msk_channel_str).centroids;
        S = segmenter(I, image_bits, obj.exp_info.params.(CHANNEL_name).(seg_type), t, mask, mask_centroids);
    else
        S = segmenter(I, image_bits, obj.exp_info.params.(CHANNEL_name).(seg_type), t);
    end
    
    % Loop over steps. params from obj.obj.exp_info. 
    seg_steps = obj.exp_info.steps.(CHANNEL_name).(seg_type);
    for i = 1:length(seg_steps)

        disp(seg_steps(i).type);
        S.(seg_steps(i).type);
    end

    %% Output. 
    % Create frame_obj and save. 
    seg_channel_name=strcat('seg_',CHANNEL_name);
    frame_obj = struct(seg_channel_name,[]);
    
    % change post-processing depending on seg_type. 
    switch seg_type
        
        case 'cells'
            
            borders = border_frame( size(S.BW) );

            % add stats to frame_obj. 
            for i = 1:length(S.stats)       
                frame_obj.(seg_channel_name).PixelIdxList{i} = S.stats(i).PixelIdxList;
                frame_obj.(seg_channel_name).centroids{i}    = S.stats(i).Centroid;

                this_cell = false(size(S.BW));
                this_cell(S.stats(i).PixelIdxList) = 1;
                frame_obj.(seg_channel_name).touches_border(i) = any(this_cell.*borders,'all');
            end

            %Add final binarized image to frame_obj for save keeping
            frame_obj.(seg_channel_name).BW = S.BW;

            %Trace boundaries.
            cBW = max(S.BW, [], 3);
            C = bwboundaries(cBW);
            if ~isempty(C)
                C = cellfun(@(x) [smooth(x(:,2)),smooth(x(:,1))],C,'uniformoutput',0);
                c_ctr = cellfun(@(x) [mean(x(:,1)),mean(x(:,2))],C,'uniformoutput',0);
                %Match centroids. 
                fobj_ctrs = cat(1,frame_obj.(seg_channel_name).centroids{:});
                D = pdist2(cat(1,c_ctr{:}),fobj_ctrs(:,1:2));
                [~,idx] = min(D);
                frame_obj.(seg_channel_name).contours=C( idx );
            end

        case 'foci'
        
            %Add final binarized image to frame_obj for save keeping
            frame_obj.(seg_channel_name).BW = S.BW;
            %Add sub-nuclear segmented nucleoli. 
            frame_obj.(seg_channel_name).foci=S.stats;
            
    end
       
    %Save frame_obj
    fname = ['frame_',sprintf('%04d',t),'.mat'];
    
    if(~exist(obj.exp_info.nuc_seg_dir,'dir'))
        mkdir(obj.exp_info.nuc_seg_dir)
    end
        
    parsave([obj.exp_info.nuc_seg_dir,fname],frame_obj,seg_channel_name)
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
function parsave(fname, frame_obj, channel_str)

if exist( fname, 'file')
    try
        D = load(fname,'frame_obj');
    catch
        D = struc();
        warning('encountered empty frame file')
    end
    
    D.frame_obj.(channel_str) = frame_obj.(channel_str);
    frame_obj = D.frame_obj;
    try
        save(fname,'frame_obj','-append');
    catch
        %For large files...
        save(fname, 'frame_obj','-v7.3','-append');
    end
else
    
    try
        save(fname,'frame_obj');
    catch
        %For large files...
        save(fname, 'frame_obj','-v7.3');
    end

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