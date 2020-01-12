% 2D Segmentation of cell nuclei. Based on: 'Object segmentation and ground
% truth in 3D embryonic imaging'. 

% This version is written for 2D images. If there's a z-stack, then just
% max-int project first. 

% 10-26-18: now takes vector of threshold values for each time frame.
function obj = segment_2D(obj, params, step, FORCE_FRAMES)


params=default_params(params);
step=default_step(step);

%% Pre-processing using bio-formats. 

%Generate reader. Use memoizer. Oif files can be really slow. Since segmentation is
%the first part of processing, there should be memo for all subsequent
%steps. 
bfInitLogging;

if step.use_specific_img    
    loci.common.DebugTools.setRootLevel('WARN');
    reader = bfGetReader();
    reader = loci.formats.Memoizer(reader,0);
    reader.setId( params.img_file );
    series = 1;    
else
    loci.common.DebugTools.setRootLevel('WARN');
    reader = bfGetReader();
    reader = loci.formats.Memoizer(reader,0);
    reader.setId( obj.exp_info.img_file );
    series = 1;
end




%Get the image size of this series. 

T = reader.getSizeT;

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

%Now run loops 
    if(parallel)
        disp('New frames to calculate:')
        disp(new_ts)
        parfor g = 1:length(new_ts)
            t = new_ts(g);
            inner_function( exp_info, Z_cell{t}, params, t, reader)
        end
    else
        disp('New frames to calculate:')
        disp(num2str(new_ts))
        for t = new_ts
            inner_function( exp_info, Z_cell{t},params,t,reader)
        end
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
function inner_function( exp_info, planes, params, t,reader)
    
%Display 
disp(['Started frame: ',num2str(t)])

%%%%%%%%%%%%%%% START PROCESSING %%%%%%%%%%%%%%%%
    image_bits     = reader.getBitsPerPixel;
    size_x = reader.getSizeX;
    size_y = reader.getSizeY;
    ZSlicesinStack = length(planes);

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
    
    
  
    % use segmenter2D to run segmentation steps. 
    S = segmenter2D(I, image_bits, params, t);
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
    
    % add stats to frame_obj. 
    for i = 1:length(S.stats)
        
        frame_obj.(channel_str).PixelIdxList{i} = S.stats(i).PixelIdxList;
        frame_obj.(channel_str).centroids{i}    = S.stats(i).Centroid;
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
