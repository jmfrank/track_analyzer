% 3D Segmentation of cell nuclei. Based on: 'Object segmentation and ground
% truth in 3D embryonic imaging'. 

% 3-09-21: Now object must contain steps/params for segmentation. 
% 7-23-18: changed the filtering step. Using diffusion in 3D code. Also,
% calculating gaussian gradient in 3D using separable filters. Now we don't
% have to loop over z. Hopefully this is a bit faster! ALSO, changing the
% thresholding step to use an absolute value rather than fraction of total
% luminance. 
% Also, keeping reader open during loop over time. Might be taking a long
% time to open and close image reader. 

function obj = segment_3D(obj, FORCE_ALL_FRAMES, debug)


%% Pre-processing using bio-formats. 

%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = bfGetReader(obj.exp_info.img_file);
series = 1; % Need to further define if multiple series in bf format. 

%Get the image size of this series. 
T = reader.getSizeT;

if nargin < 3
    debug=0;
    FORCE_ALL_FRAMES=1;
end

%% Generate im_info structure

    %Experiment info to pass along. 
    exp_info = obj.exp_info;
    % Determine which channel to perform cell segmentation on. 
    chs = fieldnames(exp_info.steps);
    CHANNEL_name=[];
    CHANNEL_id = [];
    for i = 1:length(chs)
        name=chs{i};
        if isfield(exp_info.steps.(name),'cells')
            CHANNEL_name=name;
            CHANNEL_id  =str2num(name(9:end));
            break
        end
    end
    if isempty(CHANNEL_name)
        errordlg('No cell segmentation channel found')
    end
   
        
    
%% Loop over images. Need to check if more than one image, do parallel. 



%Now check to see if any frames already exist. 
frame_files = obj.get_frame_files;

if( ~exist(frame_files{1},'file') || FORCE_ALL_FRAMES)
    new_ts = [1:T];
else
    for i = 1:length(frame_files)
        [~,fname,~] = fileparts(frame_files{i});
        t_val(i) = str2num(fname(end-3:end));
    end
    all_ts = [1:T];
    new_ts = setdiff(all_ts,t_val);
end

%If user specifies list of frames as last input, then for frames accordinly
if nargin==4
    new_ts = FORCE_FRAMES;
end


%Now run loops. Params is passed around in case things change from user. 
disp('New frames to calculate:')
disp(num2str(new_ts))
for t = new_ts
   inner_function( exp_info, CHANNEL_name, t, reader, debug);
end

%Close reader at very end. 
reader.close();

%Always clear out flags in after new segmentation. 
obj = obj.clear_flags;
obj.save; %Auto-save. 

end

%% Inner function to run 
function inner_function( exp_info, CHANNEL_name, t, reader, debug)
    
%Display 
disp(['Started frame: ',num2str(t)]);      

%%%%%%%%%%%%%%% START PROCESSING %%%%%%%%%%%%%%%%
    image_bits     = reader.getBitsPerPixel;
    CHANNEL_id     = str2num(CHANNEL_name(9:end));
    
    %Get the bio-formats image index corresponding to this z-stack:
    I = get_stack(reader, t, CHANNEL_id);
    
    %Reduce i for debugging
    if(debug)
        I = I(1:512,1:512,:);
    end

    % Check seg_dims. Max project if 2D. 
    if exp_info.params.(CHANNEL_name).cells.seg_dims == 2
        I = max(I,[],3);
    end
    
    % Mask image
    if exp_info.steps.mask_channel>0
        % Get Mask image. 
        F = obj.get_frame_files;
        load(F{1})
        msk_channel_str = ['seg_channel_',pad(num2str(exp_info.steps.mask_channel),2,'left','0')];
        mask = frame_obj.(msk_channel_str).BW;
    else
        mask = ones(size(I));
    end
    
    % Loop over steps. params from obj.exp_info. 
    S = segmenter(I, image_bits, exp_info.params.(CHANNEL_name).cells, t, mask);
    seg_steps = exp_info.steps.(CHANNEL_name).cells;
    for i = 1:length(seg_steps)

        disp(seg_steps(i).type);
        S.(seg_steps(i).type);
    end

    %% Output. 
    % Create frame_obj and save. 
    seg_channel_name=strcat('seg_',CHANNEL_name);
    frame_obj = struct(seg_channel_name,[]);
    
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
    
    %Save frame_obj
    fname = ['frame_',sprintf('%04d',t),'.mat'];
    
    if(~exist(exp_info.nuc_seg_dir,'dir'))
        mkdir(exp_info.nuc_seg_dir)
    end
        
    parsave([exp_info.nuc_seg_dir,fname],frame_obj,seg_channel_name)
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