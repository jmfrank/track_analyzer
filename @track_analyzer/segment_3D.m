% 3D Segmentation of cell nuclei. Based on: 'Object segmentation and ground
% truth in 3D embryonic imaging'. 

% Specify segmentation parameters in params, and function options in step.
% Written for single time-point images. 

% 7-23-18: changed the filtering step. Using diffusion in 3D code. Also,
% calculating gaussian gradient in 3D using separable filters. Now we don't
% have to loop over z. Hopefully this is a bit faster! ALSO, changing the
% thresholding step to use an absolute value rather than fraction of total
% luminance. 
% Also, keeping reader open during loop over time. Might be taking a long
% time to open and close image reader. 

function obj = segment_3D(obj, params, step)


% Step is optional?
if nargin < 3
    step = struct();
end

step = default_step(step);
params = default_params(params);

%% Pre-processing using bio-formats. 

%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = bfGetReader(obj.exp_info.img_file);
series = 1; % Need to further define if multiple series in bf format. 

%Get the image size of this series. 

T = reader.getSizeT;

 %Get step.debug
if(step.debug)
    params
end  

%% Generate im_info structure
    ZSlicesinStack = reader.getSizeZ;

    image_bits     = reader.getBitsPerPixel;   
    
    %Since we do parallel processing for time points, pre-compute the list
    %of plane indices for each time point. 
    Z_cell = {};
    for t = 1:T
        planes = [];
        for Z = 1:ZSlicesinStack
            planes(Z) = reader.getIndex(Z-1,params.seg_channel-1,t-1)+1;
        end
        Z_cell{t} = planes;
    end
        
    
%% Loop over images. Need to check if more than one image, do parallel. 

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
end


%Now run loops. Params is passed around in case things change from user. 
disp('New frames to calculate:')
disp(num2str(new_ts))
for t = new_ts
    params = inner_function( exp_info, Z_cell{t},params,step,t, reader);
end

%Close reader at very end. 
reader.close();

%Always clear out flags in after new segmentation. 
obj = obj.clear_flags;
obj.save;%Auto-save. 

end

%% Inner function to run 
function params = inner_function( exp_info, planes, params, step, t, reader)
    
%Display 
disp(['Started frame: ',num2str(t)])

%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%  
        

%%%%%%%%%%%%%%% START PROCESSING %%%%%%%%%%%%%%%%
    image_bits     = reader.getBitsPerPixel;
    size_x = reader.getSizeX;
    size_y = reader.getSizeY;
    ZSlicesinStack = length(planes);

    %Get image planes for this time point. 
    I = zeros(size_y,size_x,ZSlicesinStack);

    %Get the bio-formats image index corresponding to this z-stack:
    for i = 1:ZSlicesinStack
        this_plane_img = bfGetPlane(reader,planes(i));
        %If there are zeros in this plane due to weird stitching, then
        %adjust. 
        sel = this_plane_img == 0;
        if(sum(sel(:))>0)
            this_plane_img(sel) = mean( this_plane_img(~sel) );
        end
        plane_mean(i) = mean(this_plane_img(:));
        plane_std(i)  = std(double(this_plane_img(:)));
        I(:,:,i)     = this_plane_img;
    end
    
    %Reduce i for debugging
    if(step.debug)
        I = I(1:512,1:512,:);
    end
    
    S = segmenter(I, image_bits, params, t);
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
    cBW = max(S.BW, [], 3);
    C = bwboundaries(cBW);
    if ~isempty(C)
        C = cellfun(@(x) [smooth(x(:,2)),smooth(x(:,1))],C,'uniformoutput',0);
        c_ctr = cellfun(@(x) [mean(x(:,1)),mean(x(:,2))],C,'uniformoutput',0);
        %Match centroids. 
        fobj_ctrs = cat(1,frame_obj.(channel_str).centroids{:});
        D = pdist2(cat(1,c_ctr{:}),fobj_ctrs(:,1:2));
        [~,idx] = min(D);
        frame_obj.(channel_str).contours=C( idx );
    end
    
    %Save frame_obj
    fname = ['frame_',sprintf('%04d',t),'.mat'];
    
    if(~exist(exp_info.nuc_seg_dir,'dir'))
        mkdir(exp_info.nuc_seg_dir)
    end
        
    parsave([exp_info.nuc_seg_dir,fname],frame_obj,channel_str)
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

%% defaults. 
function step = default_step( step )

%List of all default parameters. 
dstep.use_specific_img=0;
dstep.FORCE_ALL_FRAMES=1;
dstep.debug=0;
dstep.subtract_bg=0;
dstep.CLAHE=0;
dstep.MeanFilter=0;
dstep.threshold_by_histogram=0;
dstep.nuc_thresh = 0;
dstep.iterative_thresholding=0;
dstep.channel=1;
dstep.merger=0;
dstep.gaussian_filter = 0;
dstep.imclose = 0;
dstep.mask = 0;
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

function params = default_params( params )



%List of all default parameters. 
dparams.view_channel=1;
dparams.channel=1;
S  = fieldnames( dparams );

for i = 1:length(S)
    
    %Check if this field exists. 
    if ~isfield(params,S{i})
        params.(S{i}) = dparams.(S{i});
        %Output this default was used. 
        disp(['Using default ',S{i},' with value: ',num2str(dparams.(S{i}))]);
    end
end



end