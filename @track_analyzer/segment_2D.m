% 2D Segmentation of cell nuclei. Based on: 'Object segmentation and ground
% truth in 3D embryonic imaging'. 

% This version is written for 2D images. If there's a z-stack, then just
% max-int project first. 

% 10-26-18: now takes vector of threshold values for each time frame.
function obj = segment_2D(obj, FORCE_ALL_FRAMES, debug)


% Step is optional?
if nargin < 3
    step = struct();
end

%params=default_params(params);
step=default_step(step);

%% Pre-processing using bio-formats. 
reader = bfGetReader(obj.exp_info.img_file);


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
if(T>1)
    parallel = 0;
else
    parallel = 0;
end

%Experiment info to pass along. 
exp_info = obj.exp_info;

%Now check to see if any frames already exist. 
frame_files = obj.get_frame_files;
if(isempty(frame_files{1}) || step.FORCE_ALL_FRAMES)
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
    
    if any( new_ts > T )
        error(['frame requested exceeds image frame count: ',num2str(T)]);
    end
end

%Now run loops. Params is passed around in case things change from user. 
disp('New frames to calculate:')
disp(num2str(new_ts))
for t = new_ts
   inner_function( exp_info, CHANNEL_name, t, reader, debug);
end

disp(['Max frames: ',num2str(T)]);

%Add params to exp_info.
%obj.exp_info.seg_params=params;
%Always clear out flags in after new segmentation. 
obj = obj.clear_flags;
%Auto-save. 
obj.save; 


end

%% Inner function to run 
function inner_function( exp_info, CHANNEL_name, params, t,reader)
    
%Display 
disp(['Started frame: ',num2str(t)])

%%%%%%%%%%%%%%% START PROCESSING %%%%%%%%%%%%%%%%
    image_bits     = reader.getBitsPerPixel;
    size_x = reader.getSizeX;
    size_y = reader.getSizeY;
    ZSlicesinStack = length(planes);
    CHANNEL_id     = str2num(CHANNEL_name(9:end));
    
    %Get image planes for this time point. 
    I = zeros(size_y,size_x,ZSlicesinStack);

    %Get the bio-formats image index corresponding to this z-stack:
    if ZSlicesinStack > 1
        if isfield(params,'z_range')
            planes = planes(params.z_range);
        end
        for i = 1:length(planes)
            this_plane_img = bfGetPlane(reader,planes(i));
            I(:,:,i)     = this_plane_img;
        end
        %Project z-dimension. 
        I = max(I,[],3);
    else
        I = bfGetPlane(reader,planes);
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
    
    % use segmenter2D to run segmentation steps. 
    S = segmenter(I, image_bits, exp_info.params.(CHANNEL_name).cells, t, mask);
    channel_str = ['channel_',pad(num2str(params.seg_channel),2,'left','0')];
    seg_steps = exp_info.steps.(channel_str).cells;
    for i = 1:length(seg_steps)
        
        disp(seg_steps(i).type);
        S.(seg_steps(i).type);
        
    end
   
    %% Output. 
    % Image channel
    channel_str = ['seg_channel_',pad(num2str(params.seg_channel),2,'left','0')];

    % Create frame_obj and save. 
    frame_obj = struct(channel_str,[]);
    
    
    borders = border_frame( size(S.BW) );

    % add stats to frame_obj. 
    for i = 1:length(S.stats)
        
        frame_obj.(channel_str).PixelIdxList{i} = S.stats(i).PixelIdxList;
        frame_obj.(channel_str).centroids{i}    = S.stats(i).Centroid;
        
        this_cell = false(size(S.BW));
        this_cell(S.stats(i).PixelIdxList) = 1;
        frame_obj.(channel_str).touches_border(i) = any(this_cell.*borders,'all');

    end
    
    
    %Add final binarized image to frame_obj for save keeping
    frame_obj.(channel_str).BW = S.BW;

    %Trace boundaries. 
    C = bwboundaries(S.BW);
    if ~isempty(C)
        C = cellfun(@(x) [smooth(x(:,2)),smooth(x(:,1))],C,'uniformoutput',0);
        c_ctr = cellfun(@(x) [mean(x(:,1)),mean(x(:,2))],C,'uniformoutput',0);
        %Match centroids. 
        D = pdist2(cat(1,c_ctr{:}),cat(1,frame_obj.(channel_str).centroids{:}));
        [~,idx] = min(D);
        frame_obj.(channel_str).contours=C( idx );
    end
        

    %% Save frame_obj as done before. 
    fname = ['frame_',sprintf('%04d',t),'.mat'];
    if(~exist(exp_info.nuc_seg_dir,'dir'))
        mkdir(exp_info.nuc_seg_dir)
    end
    parsave([exp_info.nuc_seg_dir,fname],frame_obj)
    if isempty(S.stats)
        nFound=0;
    else
        nFound=length(frame_obj.(channel_str).contours);
    end
    disp(['Finished frame:     ',num2str(t),' Found ',num2str(nFound),' cells.'])

end


%Parallel saving technique
function parsave(fname, frame_obj)
try
    save(fname,'frame_obj');
catch
    %For large files...
    save(fname, 'frame_obj','-v7.3');
end
end


%Plotting text over image
function plot_liar_text( stats, liar_list )

delete(findobj(gca,'Type','text'))
hold on
for i = 1:length(liar_list)
    ctr = stats( liar_list(i) ).Centroid;
    text(ctr(1),ctr(2),num2str( liar_list(i)),'color','g')
end
hold off

end

%% defaults. 
function step = default_step( step )

%List of all default parameters. 
dstep.use_specific_img=0;
dstep.FORCE_ALL_FRAMES=1;
S  = fieldnames( dstep );

for i = 1:length(S)
    
    %Check if this field exists. 
    if ~isfield(step,S{i})
        step.(S{i}) = dstep.(S{i});
        %Output this default was used. 
        disp(['Using default ',S{i},' with value: ',num2str(dstep.(S{i}))]);
    end
end

end
