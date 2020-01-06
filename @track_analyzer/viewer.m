% Viewing function to look at results, flagging bad cells/spots, and grouping cells together to annotate division events.
% The step structure specifies what you want to view (cells,tracks,spots...)

% This is based off imshow3D (https://www.mathworks.com/matlabcentral/fileexchange/41334-imshow3d). 

function viewer(obj, step, varargin)

%Clear out existing gui data. 
clear_app_data;

if nargin<2
    step = struct;
end

step=default_step(step);

if ~isfield(step,'AutoAdjust')
    step.AutoAdjust=0;
end

%Display tracks. 
if ~isfield(step,'tracks')
    states.tracks = 0;
end

%Display spots?
if(~isfield(step,'spots'))
    step.spots=0;
    states.spots=0;
else
    states.spots=step.spots;
end

%Parse other inputs. 
%Input parsing. 
p = inputParser;
p.addParameter('indices', [1:length(obj.tracks)]);
p.parse(varargin{:});

%% Display cell contours
if(~isfield(step,'cells'))
    step.cells = 0;
    states.cells=0;
else
    if(step.cells)
        cells = obj.getCellContours();
        states.cells=1;
    else
        states.cells=0;
    end
end


%% Default display tracks
states.group=0;
states.plot_data=0;
states.flag_tracks=0;
states.flag_cells=0;
states.flag_spots=0;
states.tracks=step.tracks;
%Save states to app data. 
setappdata(0,'states',states');

%Set drawn_cells to existing drawn_cells data. 
if isfield( obj.exp_info,'drawn_cells' )
    setappdata(0,'drawn_cells',obj.exp_info.drawn_cells )
else
    setappdata(0,'drawn_cells',[]);
end

%% Pre-processing using bio-formats.

%Figure out if data is from multiple images.  
if( iscell( obj.exp_info.img_file ))
    
    %Keep track of which img is 'on'. 
    states.img_idx = 1;
       
    %Check to make sure img_file and max_p_img are the same length. 
    if length(obj.exp_info.img_file) ~= length(obj.exp_info.max_p_img)
        error('Img file and max-p file list are different sizes.');
    end
    
    %Need to make a variable to tell us which frames are associated with which files. 
    frame2img=[];
    f_count = 0;
    n_imgs = length(obj.exp_info.max_p_img);
    Img = cell(n_imgs, 1);
    
    for e = 1:n_imgs
        
        %File name.
        fname = obj.exp_info.max_p_img{e};
        
        [reader,X,Y,Z,C,T] = bfGetReader(fname);
        
        %Keep track of frame count
        f_count = f_count + T;
        
        %This data
        frame_range = [f_count + 1: f_count + T]';
        these_frames = [1:T]';
        img_val     = e.*ones(length(frame_range),1);
        %Add to frame2img. 
        frame2img=[frame2img; frame_range,these_frames, img_val];

        this_img = bfopen(fname);
        Img{e} = zeros(Y,X,T);
        %Loop over frames, add to img cell. 
        for t = 1:T
            plane = get_planesZCT(reader,Z,step.channel,t);
            
            %If guassian filter requested. 
            if step.gaussian_filter >0
                Img{e}(:,:,t) = imgaussfilt(this_img{1}{plane,1},step.gaussian_filter);
            else
                Img{e}(:,:,t) = this_img{1}{plane,1};
            end
        end
        
        %Close reader
        reader.close();
        clear this_img

    end
    
    
else

    %Look too see if a max-p image exists. 
    if( exist( obj.exp_info.max_p_img,'file' ) )


        %Just one image. 
        fname = obj.exp_info.max_p_img;

        [reader,X,Y,Z,C,T] = bfGetReader(fname);
        frame_range = [1:T]';
        img_val = 1.*ones(length(frame_range),1);
        %Add to frame2img. 
        frame2img=[frame_range,frame_range, img_val];

        if C > 1
            %load one frame at a time. 
            Img{1} = zeros(Y,X,T);
            for t = 1:T
                plane = get_planesZCT(reader,Z,step.channel,t);
                Img{1}(:,:,t) = bfGetPlane(reader,plane);
            end
        else           
            this_img = bfopen(fname);
            Img{1} = zeros(Y,X,T);
            for t = 1:T
                plane = get_planesZCT(reader,Z,step.channel,t);
                if step.gaussian_filter > 0
                    Img{1}(:,:,t) = imgaussfilt(this_img{1}{plane,1},step.gaussian_filter);
                else
                    Img{1}(:,:,t) = this_img{1}{plane,1};
                end
            end
        end
        
        f_count = T;

        
        %Close reader
        reader.close();
        clear this_img
    %Maxp doesn't exist. 
    else
        
        %See if the img file is already a single plane. 
        [reader,X,Y,Z,C,T] = bfGetReader(obj.exp_info.img_file); 

        if Z==1
            %Load image. 
            this_img = bfopen( obj.exp_info.img_file );
            
            %Add to frame2img. 
            frame_range = [1:T]';
            img_val = 1.*ones(length(frame_range),1);
            frame2img=[frame_range,frame_range, img_val];
            
            %Compact IMG into time-sstack. Z should be 1
            Img{1} = zeros(Y,X,T);
            for t = 1:T
               planes = get_planesZCT( reader,Z,step.channel,t);
               Img{1}(:,:,t) =this_img{1}{planes,1};
            end

            %Close reader
            reader.close();
            clear this_img
        elseif T==1
            stack=get_stack(reader,T,step.channel);
            Img{1}=max(stack,[],3);
            frame2img=[1,1,1];
            Z=1
        else
            error('no maxp available');

        end
        
        f_count = T;

    end

end

% Update T
T = f_count;
%% Loading various data. 
%Add some things to appdata
setappdata(0,'track_obj',obj);

%Save dimensions to DATA
DATA.dims = [Y,X,Z,C,T];

%Get the screen size. 
screen = get(0,'ScreenSize');

%Set everything to the same figure. That way we can clear it before hand. 
MAIN = figure(7);
%MAIN.Position =  [0.75*screen(3) 0.75*screen(4) 0.24*screen(3) 0.24*screen(4)];

%Set a tag for finding this gui 
MAIN.Tag = 'Main';

clf;
set(MAIN,'color','w');

%Initialize t
t=1;
setappdata(0,'t',t);
    
%Define total number of frames. 
T = size(frame2img,1); 


%% Cell tracks and spot tracks handling.
global track_matrix track_sel_vec cell_sel_mat

%Updater function sets the above global variables. Useful for external
%control. 
tracks=[];
initialize_tracks(p.Results.indices)
       
%Initiate selection vector for spots. 
spot_sel_vec = ones(size(obj.results));

%Precalculate which spot-tracks are assigned to which tracks. 
if ~isempty( obj.spot_tracks )
    spot_track_assignment = cellfun(@(x) x(1,3), obj.spot_tracks);
end

%Look for existing flags. 
try
    
    if isfield(obj.exp_info,'flagged');
        flagged = obj.exp_info.flagged;
        disp('Found existing flags.')

        setappdata(0,'flagged',flagged)

        % Check if cells and/or tracks flagged. 
        if ~isempty(flagged.cells)

            %Loop over flagged cells. 
            for i = 1:size(flagged.cells,1)
                this_time = flagged.cells(i,1);
                this_cell_id = flagged.cells(i,2);

                % Find the track at this time. 
                this_track = track_matrix(:,this_time)==this_cell_id;

                % Set this cell off. 
                cell_sel_mat(this_track, this_time) = 0;

            end
        end
    else
        disp('no flaggs found');
    end
    
catch
    
    error('Track flags failed!')
    
end

% Cut out short tracks.
if isfield(step,'min_track_length')
    L = cellfun(@(x) size(x,1),obj.tracks);
    keep_these_tracks = find(L < step.min_track_length);
    cell_sel_mat(keep_these_tracks,:) = 0; 
end

    

%Look for existing division data
if isfield(obj.exp_info,'groups')
    if isempty(obj.exp_info.groups(1).group_id)
        groups = group_structure(1);
        group_sel_vec = zeros(size(tracks));
    else
        group_sel_vec = zeros(size(tracks));
        groups = obj.exp_info.groups;
        %Get all tracks part of groups
        for i = 1:length(groups)
            gtracks = groups(i).cell_tracks;
            group_sel_vec(gtracks) = groups(i).group_id;
        end
    end
else
    groups = group_structure(1);
    group_sel_vec = zeros(size(tracks));
end

%Colors used for on and off. 
c_vec = [1,0,0; 0,1,0];


%% Set up GUI. 
%Dimensions of image
sizes = size(Img);
M=sizes(1);
N=sizes(2);
Z=1;

global InitialCoord;

MinV = 0;
MaxV = max(Img{1}(:));
LevV = (double( MaxV) + double(MinV)) / 2;
Win = double(MaxV) - double(MinV);
WLAdjCoe = (Win + 1)/1024;
FineTuneC = [1 1/16];    % Regular/Fine-tune mode coefficients

if isa(Img{1},'uint8')
    MaxV = uint8(Inf);
    MinV = uint8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img{1},'uint16')
    MaxV = uint16(Inf);
    MinV = uint16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img{1},'uint32')
    MaxV = uint32(Inf);
    MinV = uint32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img{1},'uint64')
    MaxV = uint64(Inf);
    MinV = uint64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img{1},'int8')
    MaxV = int8(Inf);
    MinV = int8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img{1},'int16')
    MaxV = int16(Inf);
    MinV = int16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img{1},'int32')
    MaxV = int32(Inf);
    MinV = int32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img{1},'int64')
    MaxV = int64(Inf);
    MinV = int64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img{1},'logical')
    MaxV = 0;
    MinV = 1;
    LevV =0.5;
    Win = 1;
    WLAdjCoe = 0.1;
end    

SFntSz = 9;
LFntSz = 10;
WFntSz = 10;
LVFntSz = 9;
WVFntSz = 9;
BtnSz = 10;
ChBxSz = 10;

%Auto disp-range
[Rmin Rmax] = WL2R(Win, LevV);

%% Initialize image. 

figure(MAIN);
%axes('position',[0,0.2,1,0.8]), 
img_plot = subplot('Position',[0,0.2,1,0.8]);

if step.roi
    imshow(squeeze(Img{1}(y_range,x_range,t)), [Rmin Rmax]);
    %Adjust axes.
    AX = findobj( MAIN.Children, 'Type','Axes');
    AX.XLim(2) = length(x_range);
    AX.YLim(2) = length(y_range);
else
    imshow(squeeze(Img{1}(:,:,t)), [Rmin Rmax]);
    x_range = 1:X;
    y_range = 1:Y;
end

%% Deal with specified ROI. 
if step.roi
    update_roi( step.roi );
else
    shift_vec=[0,0];
end

%% Initialize image. 

if step.roi
    imshow(squeeze(Img{1}(y_range,x_range,t)), [Rmin Rmax]);
    %Adjust axes.
    AX = findobj( MAIN.Children, 'Type','Axes');
    AX.XLim(2) = length(x_range);
    AX.YLim(2) = length(y_range);
else
    imshow(squeeze(Img{1}(:,:,t)), [Rmin Rmax]);
end


%% Make text handle array. Index of entry corresponds to track id.
n_tracks= length(track_sel_vec);
all_textH = gobjects(n_tracks,1);
for i = 1:n_tracks
    all_textH(i) = text(1,1,num2str(i),'ButtonDownFcn',@textcallback,'Visible','off','UserData',i);
end

FigPos = get(MAIN,'Position');
T_pos = [50 45 uint16(FigPos(3)-100)+1 20];
time_txt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];

Wtxt_Pos = [50 20 60 20];
Wval_Pos = [110 20 60 20];
Ltxt_Pos = [175 20 45 20];
Lval_Pos = [220 20 60 20];

BtnStPnt = uint16(FigPos(3)-250)+1;
if BtnStPnt < 300
    BtnStPnt = 300;
end

Btn_Pos = [BtnStPnt 20 100 20];
ChBxFilter_Pos = [BtnStPnt-100, 20, 50, 20];
thresh_val_Pos  = [BtnStPnt-50, 20, 50, 20];
thresh_text_Pos = [BtnStPnt 20 45 20];

%Time label. 
time_txthand = uicontrol('Style', 'text','Position', time_txt_Pos,'String',sprintf('Time# %d / %d',t, T), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
ltxthand = uicontrol('Style', 'text','Position', Ltxt_Pos,'String','Level: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
wtxthand = uicontrol('Style', 'text','Position', Wtxt_Pos,'String','Window: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', WFntSz);
lvalhand = uicontrol('Style', 'edit','Position', Lval_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @WinLevChanged);
wvalhand = uicontrol('Style', 'edit','Position', Wval_Pos,'String',sprintf('%6.0f',Win), 'BackgroundColor', [1 1 1], 'FontSize', WVFntSz,'Callback', @WinLevChanged);
if(T>1)
    Thand    = uicontrol('Style', 'slider','Min',1,'Max',T,'Value',t,'SliderStep',[1/(T-1) 10/(T-1)],'Position', T_pos,'Callback', @TimeSlider);
else
    Thand    = uicontrol('Style','text','string','Single time-point','Position',T_pos);
end

%gui objects for thresholding spots
ChBxFilter_hand = uicontrol('Style','checkbox','Position',ChBxFilter_Pos,'String','Filter','BackgroundColor',[0.8 0.8 0.8], 'FontSize', ChBxSz);
thresh_text_hand = uicontrol('Style', 'text','Position', thresh_text_Pos,'String','Threshold: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
thresh_val_hand = uicontrol('Style', 'edit','Position', thresh_val_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback',@DrawSpots);

%Add slider to guidata. 
gdata.slider=Thand;
gdata.reset_flags = @reset_flags;
gdata.cells_callback = @cells_callback;
gdata.update_roi     = @update_roi;
gdata.initialize_tracks  = @initialize_tracks;
guidata(MAIN, gdata)

set(MAIN, 'WindowScrollWheelFcn', @mouseScroll);
set(MAIN, 'ButtonDownFcn', @mouseClick);
%set(get(gca,'Children'),'ButtonDownFcn', @mouseClick); << not sure what
%this was doing ? 
set(MAIN,'WindowButtonUpFcn', @mouseRelease)
set(MAIN,'ResizeFcn', @figureResized)
set(MAIN,'keypressfcn',@fh_kpfcn);
%Run track/spot drawing. Needs to be after buttondown Fcn set on Main
%because text are children of Main? 
DrawTracks
DrawCells(states)
DrawSpots


%% tool box
%Need to tell gui about existing groups. 
setappdata(0,'groups',groups);
setappdata(0,'group_sel_vec',group_sel_vec);

%Delete existing tool box. 
delete( findobj('Tag','TOOL_BOX') );
TOOL = TOOL_BOX( );
TOOL.Visible='on';

%% guidata functions - created to allow external user control. 

% -=< Reset flags matrix >=-
    function reset_flags()
        
        %Cell selection matrix. 1 means active, 0s means inactive. 
        cell_sel_mat = false(size(track_matrix));
        cell_sel_mat( track_matrix > 0 ) = 1;
        
        flagged = getappdata(0,'flagged');
        states  = getappdata(0,'states');
        
        if states.flag_cells || states.flag_tracks
            
            flagged.cells=[];
        
        elseif states.flag_spots
            
            flagged.spots=[];
        end
        
        setappdata(0,'flagged',flagged);
        DrawTracks
        DrawCells(states)
        DrawSpots  
    end


% If user decides they really want to look at cells and they weren't
% loaded. 
    function cells_callback()
        cells = obj.getCellContours();
    end

% Initialize the cell track variables. 
    function initialize_tracks( indices )
        
        %Default is to turn all tracks 'on'. 
        if nargin < 0
            indices = [1:length(obj.tracks)];
        end
       
        tracks=obj.tracks;
        
        %To save time in plotting, make a reference matrix to tell us which tracks
        %are part of which frame. While looping, also fill out cell_sel data. 
        track_matrix = zeros(length(tracks),T);
        
        for i = indices
            ts = tracks{i}(:,1);
            cell_ids = tracks{i}(:,2);
            track_matrix(i,ts) = cell_ids';
        end

        %Initiate selection vector for cell tracks. 
        track_sel_vec = ones(size(tracks));

        %Cell selection matrix. 1 means active, 0s means inactive. 
        cell_sel_mat = false(size(track_matrix));
        cell_sel_mat( track_matrix > 0 ) = 1;      

    end

% Allow adjusting the current selection of visible tracks. 

    function update_track_selection( indices )
        
        off_sel = setdiff(1:length(obj.tracks),indices);
        track_matrix = track_matrix_og;
        
        track_matrix( off_sel, : ) = 0;
        
    end

% Update frame. 
    function update_frame
        
        %Get image index. 
        img_idx = frame2img(t,3);
        frame_idx = frame2img(t,2);
        
        H=findobj(MAIN,'type','Image');
        
        if step.roi
            set(H,'cdata',squeeze(Img{img_idx}(y_range,x_range,frame_idx)));
            %Adjust axes.
            AX.XLim(2) = length(x_range);
            AX.YLim(2) = length(y_range);
        else
            set(H,'cdata',squeeze(Img{img_idx}(:,:,frame_idx)));
        end
        

    end
        
% Update ROI. input roi must be [lower_x, lower_y, width, height]
    function update_roi( new_roi )
        
        %If no input variable, reset to no ROI. 
        if nargin < 1
            step.roi = 0;
            shift_vec = [0,0];
            AX = findobj( MAIN.Children, 'Type','Axes');
            AX.XLim(2) = X;
            AX.YLim(2) = Y;
        else            

            x_range = new_roi(1): new_roi(1) + new_roi(3) - 1;
            y_range = new_roi(2): new_roi(2) + new_roi(4) - 1;
            %Check size of image. 
            sel_x = x_range > 0 & x_range < X;
            sel_y = y_range > 0 & y_range < Y;
            
            if sum(sel_x)==0 | sum(sel_y)==0
                disp('bad roi')
                return
            end
            x_range = x_range(sel_x);
            y_range = y_range(sel_y);

            step.roi = new_roi;

            AX = findobj( MAIN.Children, 'Type','Axes');
            AX.XLim(2) = length(x_range);
            AX.YLim(2) = length(y_range);
            
            %%% There might be a bug when the plot  was zoomed in and ROI was
            %%% updated. 



            %Shift if using ROI. 
            shift_vec = [-step.roi(1),-step.roi(2)];
        end
        
    end

%% Original GUI functions from imshow3D. 
% -=< Figure resize callback function >=-
    function figureResized(object, eventdata)
        FigPos = get(MAIN,'Position');
        T_Pos        = [50 45 uint16(FigPos(3)-100)+1 20];
        time_txt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
        BtnStPnt = uint16(FigPos(3)-100)+1;
        if BtnStPnt < 300
            BtnStPnt = 300;
        end
        Btn_Pos = [BtnStPnt 20 100 20];
        
        ChBxFilter_Pos = [BtnStPnt-100, 20, 50, 20];
        thresh_val_Pos  = [BtnStPnt-50, 20, 50, 20];
        thresh_text_Pos = [BtnStPnt 20 45 20];
        
        set(time_txthand,'Position',time_txt_Pos);
        set(Thand,'Position', T_Pos);
        %set(Btnhand,'Position', Btn_Pos);
        
        set(ltxthand,'Position', Ltxt_Pos);
        set(wtxthand,'Position', Wtxt_Pos);
        set(lvalhand,'Position', Lval_Pos);
        set(wvalhand,'Position', Wval_Pos);

        ChBxFilter_Pos = [BtnStPnt-100, 20, 50, 20];
        thresh_val_Pos  = [BtnStPnt-50, 20, 50, 20];
        set(thresh_text_hand,'Position',thresh_text_Pos);
        set(ChBxFilter_hand,'Position',ChBxFilter_Pos);
        set(thresh_val_hand,'Position',thresh_val_Pos);
    end

% -=< Time slider callback function >=-
    function TimeSlider(object,event)
        t = round(get(object,'Value'));
        setappdata(0,'t',t);
        
        update_frame
        
        if step.AutoAdjust
            AutoAdjust
        end

        if T > 1
            set(time_txthand, 'String', sprintf('Time# %d / %d',t,T));
        else
            set(time_txthand, 'String', '2D image');
        end
        DrawTracks
        %Add drawspots here so its' called everytime DrawTracks is. 
        DrawCells(states)
        DrawSpots
    end

% -=< Mouse scroll wheel callback function >=-
    function mouseScroll (object, eventdata)
        MAIN;
        UPDN = eventdata.VerticalScrollCount;
        t = t - UPDN;
        if (t < 1)
            t = 1;
        elseif (t > T)
            t = T;
        end
        setappdata(0,'t',t);

        if T > 1
            set(time_txthand, 'String', sprintf('Time# %d / %d',t,T));
        else
            set(time_txthand, 'String', '2D image');
        end
        
        update_frame
        
        set(Thand,'Value',t);
        DrawTracks
        
        %Add drawspots here so its' called everytime DrawTracks is. 
        DrawCells(states)
        DrawSpots
    end

% -=< Mouse button released callback function >=-
    function mouseRelease (object,eventdata)
        set(MAIN, 'WindowButtonMotionFcn', '')
    end

% -=< Mouse click callback function >=-
    function mouseClick (object, eventdata)
        MouseStat = get(gcbf, 'SelectionType');
        if (MouseStat(1) == 'a')        %   RIGHT CLICK
            InitialCoord = get(0,'PointerLocation');
            set(MAIN, 'WindowButtonMotionFcn', @WinLevAdj);
        end
    end

% -=< Window and level mouse adjustment >=-
    function WinLevAdj(varargin)
        PosDiff = get(0,'PointerLocation') - InitialCoord;
        %Win = Win + PosDiff(1) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        %LevV = LevV - PosDiff(2) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        
        if (Win < 1)
            Win = 1;
        end
        
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
        InitialCoord = get(0,'PointerLocation');
    end

% -=< Window and level text adjustment >=-
    function WinLevChanged(varargin)

        LevV = str2double(get(lvalhand, 'string'));
        Win = str2double(get(wvalhand, 'string'));
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
    end

% -=< Window and level to range conversion >=-
    function [Rmn Rmx] = WL2R(W,L)
        Rmn = L - (W/2);
        Rmx = L + (W/2);
        if (Rmn >= Rmx)
            Rmx = Rmn + 1;
        end
    end
        
% -=< Window and level auto adjustment callback function >=-
    function AutoAdjust(object,eventdata)
        %Image handle. 
        h= findobj(MAIN,'Type','Image');
        dat = h.CData(:);
        Win = double(max(dat)-min(dat));
        Win (Win < 1) = 1;
        LevV = double(min(dat) + (Win/2));
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
    end


%% Drawing detected spot centroids on top of frame
    function DrawTracks(object, eventdata)
        %Update states. 
        states = getappdata(0,'states');
        MAIN;
        %Check if tracks exist. 
        if isempty(tracks)
            return
        end
        
        if ~states.tracks
            %Check if there are some text objects. 
            [all_textH.Visible]=deal('off');
            return
        end
        
        MAIN;
        
        %Use track matrix and cell selection matrix to figure out which
        %tracks are in this frame.
        on_selection = track_matrix(:,t) > 0 & cell_sel_mat(:,t);
        on_selection_flagged = track_matrix(:,t) > 0 & ~cell_sel_mat(:,t);
        off_selection = ~(track_matrix(:,t) > 0); 
        
        %Update positions for any object on this frame. 
        POS = cellfun(@(x) x(find(x(:,1)==t),3:4)+shift_vec,tracks,'UniformOutput',0);
        [all_textH(find( track_matrix(:,t) > 0 )).Position] = POS{on_selection | on_selection_flagged };
        
        %Update colors. 
        [all_textH(  on_selection ).Color]=deal(c_vec(2,:));
        [all_textH(  on_selection_flagged ).Color]=deal(c_vec(1,:));
        
        %set visibility. 
        [all_textH( on_selection | on_selection_flagged ).Visible] = deal('on');
        [all_textH( off_selection ).Visible] = deal('off' );
       
        
        %Grouping data
        if states.group
            group_sel_vec = getappdata(0,'group_sel_vec');
            U = unique( group_sel_vec );
            hlist = findobj(TOOL,'Tag','GroupList');
            list_string = hlist.String;
            if iscell(list_string)
                n_colors = length(hlist.String);
            else
                n_colors = 1;
            end
           
            group_colors = linspecer(n_colors);
            
            div_tracks = find( group_sel_vec > 0 & track_matrix(:,t)' > 0 );                            
            for i = div_tracks
                %c = find(U==group_sel_vec(i))-1;
                %all_textH(i).Color=group_colors(c,:);
                all_textH(i).Color= group_colors( group_sel_vec(i),:);
            end
        end
        
        hold off
        
        %Plotting data. 
        if states.plot_data
            %Update the indicator
            h = getappdata(0,'plot_time_indicator');
            if isempty(h)
                return
            end
            plot_signal = getappdata(0,'plot_signal');
            set(h(1),'XData',[t,t]);
            set(h(2),'XData',t);
            t_idx = plot_signal(:,1)==t;
            set(h(2),'YData',plot_signal(t_idx,2));
            %Back to main
            figure(MAIN);
        end
    end

%% Text click call back function. Define so we can use global sel_vec...
    function textcallback(hobj, event)
        MAIN;
        %Get current gui states.
        states=getappdata(0,'states');
        
        %Get current flagged states. 
        flagged = getappdata(0,'flagged');
        if isempty(flagged)
            flagged=struct('cells',[]);
        end

        if ~isfield(flagged,'cells')
            flagged.cells=[];
        end
            
        %Get track id from userdata.
        track_id = hobj.UserData;
        
        %Find cell id for the frame. 
        cell_id  = track_matrix(track_id,t);
        
        %Get all time points for this track. 
        t_sel = track_matrix(track_id,:) > 0;
        
        %If we are flagging tracks, change all cells in track
        if states.flag_tracks 
            
            cell_status = cell_sel_mat(track_id,t);
            
            %Switch status. 
            if cell_status 
                %Switch off
                cell_sel_mat(track_id,t_sel) = 0; 
                hobj.Color = c_vec(1,:);
                track_sel_vec(track_id) = 0;
                these_frames = find(t_sel);
                these_cell_ids = track_matrix(track_id,t_sel);
                
                %Add cells to flagged structure. 
                flagged.cells = [flagged.cells; these_frames', these_cell_ids']; 
            else
                %Switch on. 
                cell_sel_mat(track_id, t_sel) = 1;
                hobj.Color = c_vec(2,:);
                track_sel_vec(track_id) = 1;
                these_frames = find(t_sel);
                these_cell_ids = track_matrix(track_id, t_sel);
                
                %Remove from flagged structure. 
                flagged.cells = setdiff(flagged.cells,[these_frames',these_cell_ids'],'rows');
            end

            %Note the time. 
            flagged.DateString = datestr(clock);
            %Add flagged data to appdata
            setappdata(0,'flagged',flagged);
            
        elseif states.flag_cells
            
            cell_status = cell_sel_mat(track_id,t);
            
            %Switch status. 
            if cell_status 
                %Switch off
                cell_sel_mat(track_id,t) = 0; 
                hobj.Color = c_vec(1,:);
                flagged.cells = [flagged.cells; t, track_matrix(track_id,t)];
            else
                %Switch on
                cell_sel_mat(track_id,t) = 1;
                hobj.Color = c_vec(2,:);
                col_k = flagged.cells(:,1)==t & flagged.cells(:,2) == track_matrix(track_id,t);
                flagged.cells = flagged.cells(~col_k,:);
            end

            %Note the time. 
            flagged.DateString = datestr(clock);
            %Add flagged data to appdata
            setappdata(0,'flagged',flagged);
            
        elseif states.group
            
            group_sel_vec = getappdata(0,'group_sel_vec');
            
            %Ignore if the current track is flagged. 
            if ~track_sel_vec(track_id) 
                return
            end
            
            %This has yet to be assigned. Get the current division
            %track id selected in the toolbox gui.            
            curr_group_id = getappdata(0,'group_id');
            
            if isempty(curr_group_id)
                disp('Add group!')
                return
            end
            
            %Get current group data
            G = getappdata(0,'groups');
            all_ids = [G.group_id];
            this_entry = find(all_ids==curr_group_id);
            these_group_tracks = G(this_entry).cell_tracks;
            
            %Check the division_selection vector. If the track was already
            %selected for this group, switch to zero. 
            if group_sel_vec(track_id) == curr_group_id
                %Setting to zero local vector
                group_sel_vec(track_id) = 0;
                %Changing color
                hobj.Color = c_vec(2,:);
                
                %Need to remove this track from this group. 
                sel = these_group_tracks ~= track_id;
                G(this_entry).cell_tracks = these_group_tracks(sel);
                
            %If the track was unselected, then add to group data. 
            elseif  group_sel_vec(track_id) == 0
                
                %Set this division selection vector
                group_sel_vec(track_id) = curr_group_id;
               
                %Figure out coloring. 
                U = unique( group_sel_vec );
                hlist = findobj(TOOL,'Tag','GroupList');   
                list_string = hlist.String;
                if iscell(list_string)
                    n_colors = length(hlist.String);
                else
                    n_colors = 1;
                end
                
                group_colors = linspecer(n_colors);

                %Now change color
                c= find(U==curr_group_id)-1;
                hobj.Color = group_colors(c,:);
                
                %Append to group data
                G(this_entry).cell_tracks = [these_group_tracks, track_id];
                
            end
            
            %Update group data
            setappdata(0,'groups',G);
            setappdata(0,'group_sel_vec',group_sel_vec);
        else
                       
            % DO NOTHING. 
            
        end
                
        %If plot state is on, plot data associated with this track. 
        if states.plot_data

            data_plot = newFigure(8);
            hold on
            %Plot location is wrt MAIN
            %main_pos = MAIN.Position;
            %data_plot.Position =[main_pos(1) main_pos(2)-0.3*screen(4) 0.2*screen(3) 0.2*screen(4)];
        
            
            this_track = tracks{ track_id };
            data =obj.get_track_data( 'track', this_track );
            sig = cat(1,data.nuc_mean); %./ cat(1,data.cyto_mean);
            %sig = [data.nuc_area]';
            title(['Track: ',num2str(track_id)]);
            %Plot the trace. 
            plot(this_track(:,1),sig,'-*','linewidth',4,'color','k');
            
            %Plot a marker / line for the current time. 
            curr_idx = this_track(:,1)==t;

            box off
            %ylim([0.5,2])
            yvals=ylim;
            line_color=[148,47,142]./255;
            plot_time_indicator(1) = plot([t,t],yvals,'-','linewidth',5,'color',[line_color,0.6]);
            plot_time_indicator(2) = plot(t,sig(curr_idx),'.r','markersize',20);
            
            hold off
            setappdata(0,'plot_time_indicator',plot_time_indicator);
            setappdata(0,'plot_signal',[this_track(:,1),sig]);
            
            %Back to main
            figure(MAIN);
        end
                
    end

%Spot circle call back function. Upon click... 
    function circ_callback(obj,event)
    
        states = getappdata(0,'states');
        
        if states.flag_spots
            
            fit_id = obj.UserData;

            %Switch status. 
            if( spot_sel_vec( fit_id ) )
                %Switch off
                spot_sel_vec( fit_id ) = 0;
                obj.Color = c_vec(1,:);
            else
                %Switch on
                spot_sel_vec( fit_id ) = 1;
                obj.Color = c_vec(2,:);
            end
        end

        
    end

%% Drawing spots
    function DrawSpots
        %Remove old spots. 
        delete( findobj(gca,'Type','hggroup') );
        
        if ~states.spots
            return
        end
        
        MAIN;

        if( get(ChBxFilter_hand,'Value'))
            %Get filter value. 
            thresh_val =str2double(get(thresh_val_hand,'string'));
        else
            thresh_val = 0;
        end
        
        %Get cell tracks that are currently active. 
        on_tracks = find( track_matrix(:,t) > 0 & cell_sel_mat(:,t) )';   
        
        %Looping over track IDs
        for i = on_tracks
            
            %Look for spot_tracks assigned to this cell track. 
            s_tracks_ids = find( spot_track_assignment == i );
            
            if(isempty(s_tracks_ids))
                continue
            end
            
            %Figure out if any spot tracks exist at time t for cell track i. 
            these_tracks = obj.spot_tracks( s_tracks_ids );
            
            %Most likely just one track... but loop in case. 
            for j = 1:length(these_tracks)
                
                these_times = these_tracks{j}(:,1);
                id = find(these_times==t);
                
                if(~isempty(id))
                    %Get exact result id using spot track of this cell track. 
                    result_id = these_tracks{j}(id,2);
                    int = cat(1,obj.results(result_id).sum_int);
                    pos = cat(1,obj.results(result_id).pos);

                    %Check if intensity above threshold. 
                    if( int < thresh_val )
                        continue
                    end
                    
                    %Plot as distinct object.
                    if( spot_sel_vec( result_id ))
                        SC=[125,42,142]/255; %Purple...
                        SC=[172,156,255]/255;
                        SC=[113,191,110]./255;

                        h = viscircles(pos([2,1])+shift_vec, 5, 'color', SC, 'EnhanceVisibility',0,'linewidth',4);
                    else
                        h = viscircles(pos([2,1])+shift_vec, 5, 'color', 'r');
                    end
                    %Setting some properties of the circle so it's
                    %clickable
                    h.Children(1).UserData = result_id;
                    h.Children(1).ButtonDownFcn = @circ_callback;
                end
            end
        end
        
    end

%% Drawing cell contours. 
    function DrawCells( states )
        %Remove contours. 
        delete( findobj(MAIN,'type','line') ) 
        
        if(~states.cells)
            return
        end
            
        if ~exist('cells','var')
            cells_callback;
        end     
        
        MAIN;
        %Figure out dimensions. 
        if DATA.dims(3) == 1
            dim = 2; 
        else
            dim=3;
        end
        
        %Figure out which cells to draw. 
        if isfield(step,'selected_tracks')
            track_ids=step.selected_tracks;
        else
            track_ids = [1:length(tracks)];
        end
        
        %Find union of selected tracks and track_matrix data
        on_tracks = find( track_matrix(:,t) > 0 & cell_sel_mat(:,t) )';
        on_tracks = intersect(on_tracks,track_ids);        
        
        hold on
        if dim==2
            %Find the index in cells. 
            for p_ = 1:length(on_tracks)
                
                frames = tracks{on_tracks(p_)}(:,1);
                cell_id = tracks{on_tracks(p_)}(frames==t, 2);
                
                %Plot contour. 
                line_color=[148,47,142]./255;
                %line_color=[113,191,110]./255;
                %line_color=[0,0,0];
                plot(cells{t}{cell_id}(:,1)+shift_vec(1),cells{t}{cell_id}(:,2)+shift_vec(2),':','linewidth',6,'color',[line_color,0.9]);
            end
            
        elseif dim==3
            
            %Project and make convhull? 
            disp('missing code')
        end
        
        
        %If there are drawn cells, then plot ROIs. 
        drawn_cells = getappdata(0,'drawn_cells');
        if ~isempty( drawn_cells )
            hold on
            
            %Check if any drawn_cells are for this frame. 
            drawn_frames =  [drawn_cells.t];
            idx = find( drawn_frames == t );
            for c = 1:length(idx)
                
                pos = drawn_cells( idx(c) ).vector;
                plot(pos(:,1),pos(:,2),'linewidth',2,'color','w')
            end
            hold off
        end
        
        
        hold off

    end

% -=< Arrow key changes thresh value  >=-

    function [] = fh_kpfcn(H,E)          
    % Figure keypressfcn
    switch E.Key
        case 'd'
            
            TB_gd = guidata(TOOL);            
            TB_gd.draw_cells()
            
        
        case 'p'
            
            pan;
            
            
        %otherwise  
    end
    end

end

%% Remove all app data. 
function clear_app_data()

D = getappdata(0);

fields = fieldnames(D);

    %Loop
    for d = 1:length(fields)

        rmappdata(0,fields{d});
    end

end

%% Subfunction for getting the requested image.
function Img = frame_loader(reader, t, channel_id)

%Loop over z-stack. 
Z= reader.getSizeZ;
X= reader.getSizeX;
Y= reader.getSizeY;
Img = zeros(Y,X,Z);

    %Loop over z only
    for z = 1:Z
        plane_id = reader.getIndex(z-1,channel_id-1,t-1)+1;
        Img(:,:,z) = bfGetPlane(reader,plane_id);
    end

end

%% Key pressing. 

