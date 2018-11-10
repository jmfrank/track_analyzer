% Viewing program to quickly look at results. Uses max-int projection of
% raw image and plots Track id and detected nascent spots. 

% This is based off imshow3D. 
function obj = results_viewer_V2(obj, params,step,tracks)

if(nargin<2)
    params.view_channel=1;
end
if(~isfield(params,'view_channel'))
    params.view_channel = 1;
end


%Display spots?
if(~isfield(step,'spots'))
    step.spots=0;
    states.spots=0;
end


%Display cell contours
if(~isfield(step,'cells'))
    step.cells = 0;
    stats.cells=0;
    
else
    if(step.cells)
        %Load frame obj. 
        seg_files = obj.get_frame_files;
        %Loop and load. 
        disp('loading frame objs')
        for i=1:length(seg_files)
            DATA.seg(i) = load( seg_files{i},'frame_obj');
        end
        states.cells=1;
    else
        states.cells=0;
end

%Default display tracks
states.tracks=1;

%Save states to app data. 
setappdata(0,'states',states');

%% Pre-processing using bio-formats.

%Look too see if a max-p image exists. 
if( isfield( obj.exp_info,'max_p_img' ) )
   
    %Copy to local temp
    %[a,b,ex] = fileparts( obj.exp_info.max_p_img);
    %temp_loc = ['/home/jan/TEMP/',b,ex]
    %copyfile(obj.exp_info.max_p_img);
    
    %Load image. 
    IMG = bfopen( obj.exp_info.max_p_img );

    %Get the reader. 
    [reader,X,Y,Z,C,T] = bfGetReader(obj.exp_info.max_p_img);
    
    Img = zeros(Y,X,T);
    %Compact IMG into time-stack. Z should be 1....
    for t = 1:T
       planes = get_planes( reader,Z,params.view_channel,t);
       Img(:,:,t) =IMG{1}{planes,1};
    end
    
else
    
   %See if the img file is already a single plane. 
   [reader,X,Y,Z,C,T] = bfGetReader(obj.exp_info.img_file);
   
   if(Z==1)
        %Load image. 
        IMG = bfopen( obj.exp_info.img_file );

        Img = zeros(Y,X,T);
        %Compact IMG into time-stack. Z should be 1....
        for t = 1:T
           planes = get_planes( reader,Z,params.view_channel,t);
           Img(:,:,t) =IMG{1}{planes,1};
        end
   else
       %Go ahead and compute downsampled maxp
        params.seg_channel = params.view_channel;
        obj = obj.create_ds_time_series(params);
        obj.save();  

        %Load image. 
        IMG = bfopen( obj.exp_info.max_p_img );

        %Get the reader. 
        [reader,X,Y,Z,C,T] = bfGetReader(obj.exp_info.max_p_img);

        Img = zeros(Y,X,T);
        %Compact IMG into time-stack. Z should be 1....
        for t = 1:T
           planes = get_planes( reader,Z,params.view_channel,t);
           Img(:,:,t) =IMG{1}{planes,1};
        end

   end
end    

%Close reader
reader.close();

%Add some things to appdata
setappdata(0,'track_obj',obj);


%Save dimensions to DATA
DATA.dims = [Y,X,Z,C,T];

%Get the screen size. 
screen = get(0,'ScreenSize');


%Set everything to the same figure. That way we can clear it before hand. 
MAIN = figure(7);
MAIN.Position =  [0.75*screen(3) 0.75*screen(4) 0.24*screen(3) 0.24*screen(4)];

%Set a tag for finding this gui 
MAIN.Tag = 'Main';

clf;
set(MAIN,'color','w');
%Initialize t
t=1;
if(nargin < 3)
    channel_id = 1;
end
    
%Define some variables
T = size(Img,3); 

%Tracks_curr is currently used tracks. 
if nargin >3
    tracks_curr =tracks;
else
    tracks_curr = obj.tracks;
end

%To save time in plotting, make a reference matrix to tell us which tracks
%are part of which frame. 
track_matrix = zeros(length(tracks_curr),T);
track_starts = zeros(length(tracks_curr),1);
for i = 1:length(tracks_curr)
    
    ts = tracks_curr{i}(:,1);
    track_matrix(i,ts) = 1;
    track_starts(i) = ts(1);
end


%Initiate selection vector for cell tracks. 
sel_vec = ones(size(tracks_curr));
%Initiate selection vector for spots. 
spot_sel_vec = ones(size(obj.results));

%Look for existing flags. 
try
    flags = obj.exp_info.flagged.tracks;
    %Go ahead and flag these tracks. 
    sel_vec(flags) = 0;
end

%Look for existing division data
if isfield(obj.exp_info,'groups')
    if isempty(obj.exp_info.groups(1).group_id)
        groups = group_structure(1);
        group_sel_vec = zeros(size(tracks_curr));
    else
        group_sel_vec = zeros(size(tracks_curr));
        groups = obj.exp_info.groups;
        %Get all tracks part of groups
        for i = 1:length(groups)
            tracks = groups(i).cell_tracks;
            group_sel_vec(tracks) = groups(i).group_id;
        end
    end
else
    groups = group_structure(1);
    group_sel_vec = zeros(size(tracks_curr));
end


%Colors used for on and off. 
c_vec = [1,0,0; 0,1,0];

%% TOOL box
%Need to tell gui about existing groups. 
setappdata(0,'groups',groups);
setappdata(0,'group_sel_vec',group_sel_vec);

TOOL = TOOL_BOX( );


%% Set up. 
%Dimensions of image
sizes = size(Img);
M=sizes(1);
N=sizes(2);
Z=1;

global InitialCoord;

MinV = 0;
MaxV = max(Img(:));
LevV = (double( MaxV) + double(MinV)) / 2;
Win = double(MaxV) - double(MinV);
WLAdjCoe = (Win + 1)/1024;
FineTuneC = [1 1/16];    % Regular/Fine-tune mode coefficients

if isa(Img,'uint8')
    MaxV = uint8(Inf);
    MinV = uint8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint16')
    MaxV = uint16(Inf);
    MinV = uint16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint32')
    MaxV = uint32(Inf);
    MinV = uint32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint64')
    MaxV = uint64(Inf);
    MinV = uint64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int8')
    MaxV = int8(Inf);
    MinV = int8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int16')
    MaxV = int16(Inf);
    MinV = int16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int32')
    MaxV = int32(Inf);
    MinV = int32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int64')
    MaxV = int64(Inf);
    MinV = int64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'logical')
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
figure(MAIN);
%axes('position',[0,0.2,1,0.8]), 
img_plot = subplot('Position',[0,0.2,1,0.8]);
imshow(squeeze(Img(:,:,t,:)), [Rmin Rmax])


FigPos = get(MAIN,'Position');
T_pos = [50 45 uint16(FigPos(3)-100)+1 20];
time_txt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];

Wtxt_Pos = [50 20 60 20];
Wval_Pos = [110 20 60 20];
Ltxt_Pos = [175 20 45 20];
Lval_Pos = [220 20 60 20];
Track_Pos =[280 20 60 20];

BtnStPnt = uint16(FigPos(3)-250)+1;
if BtnStPnt < 300
    BtnStPnt = 300;
end
Btn_Pos = [BtnStPnt 20 100 20];
DoneBtn_Pos = [BtnStPnt+110 20 100 20];
ChBxFilter_Pos = [BtnStPnt-100, 20, 50, 20];
thresh_val_Pos  = [BtnStPnt-50, 20, 50, 20];
thresh_text_Pos = [280 20 45 20];

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

Btnhand = uicontrol('Style', 'pushbutton','Position', Btn_Pos,'String','????', 'FontSize', BtnSz, 'Callback' , @LinkTracks);
DoneBtnhand = uicontrol('Style', 'pushbutton','Position', DoneBtn_Pos,'String','Finished', 'FontSize', BtnSz);

%gui objects for thresholding spots
ChBxFilter_hand = uicontrol('Style','checkbox','Position',ChBxFilter_Pos,'String','Filter','BackgroundColor',[0.8 0.8 0.8], 'FontSize', ChBxSz);
thresh_text_hand = uicontrol('Style', 'text','Position', thresh_text_Pos,'String','Threshold: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
thresh_val_hand = uicontrol('Style', 'edit','Position', thresh_val_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback',@DrawSpots);

%Precalculate which spot-tracks are assigned to which cells. 
if(step.spots)
    fit_ids = cellfun(@(x) x(1,2), obj.spot_tracks);
    spot_track_assignment = cat(1,obj.results(fit_ids).cell_id);
end

%Run track/spot drawing. 
DrawTracks

set (MAIN, 'WindowScrollWheelFcn', @mouseScroll);
set (MAIN, 'ButtonDownFcn', @mouseClick);
set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
set(MAIN,'WindowButtonUpFcn', @mouseRelease)
set(MAIN,'ResizeFcn', @figureResized)


% -=< Figure resize callback function >=-
    function figureResized(object, eventdata)
        FigPos = get(MAIN,'Position');
        T_Pos        = [50 45 uint16(FigPos(3)-100)+1 20];
        time_txt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
        BtnStPnt = uint16(FigPos(3)-250)+1;
        if BtnStPnt < 300
            BtnStPnt = 300;
        end
        Btn_Pos = [BtnStPnt 20 100 20];
        DoneBtn_Pos = [BtnStPnt + 110 20 100 20];
        set(time_txthand,'Position',time_txt_Pos);
        set(Thand,'Position', T_Pos);
        set(Btnhand,'Position', Btn_Pos);
        
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

% -=< Slice slider callback function >=-
    function SliceSlider (hObj,event)
        S = round(get(hObj,'Value'));
        H = findobj(img_plot,'Type','Image');
        set(H,'cdata',squeeze(Img(:,:,S)));
        %set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
        caxis([Rmin Rmax])

        Cobj = findobj(img_plot,'Type','contour');
        delete(Cobj);
        BW_plane = BW(:,:,S);
        subplot(img_plot)
        hold on
        contour(BW_plane,[0.5,0.5],'color','w');
        hold off
    end

% -=< Time slider callback function >=-
    function TimeSlider(object,event)
        t = round(get(object,'Value'));
        setappdata(0,'t',t);
        H = findobj(img_plot,'Type','Image');
        set(H,'cdata',squeeze(Img(:,:,t,:)))
        if T > 1
            set(time_txthand, 'String', sprintf('Time# %d / %d',t,T));
        else
            set(time_txthand, 'String', '2D image');
        end
        DrawTracks
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
        if T > 1
            set(time_txthand, 'String', sprintf('Time# %d / %d',t,T));
        else
            set(time_txthand, 'String', '2D image');
        end
        H = findobj(img_plot,'Type','Image');
        set(Thand,'Value',t);
        set(H,'cdata',squeeze(Img(:,:,t,:)))
        DrawTracks
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
        Win = double(max(Img(:))-min(Img(:)));
        Win (Win < 1) = 1;
        LevV = double(min(Img(:)) + (Win/2));
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
    end
        

%Drawing detected spot centroids on top of frame
    function DrawTracks(object, eventdata)
        MAIN;
        %Update states. 
        states = getappdata(0,'states');
        %Check status 
        FLAG_STATE = getappdata(0,'FLAG_STATE');
        CONN_STATE = getappdata(0,'CONN_STATE');
        PLOT_STATE = getappdata(0,'PLOT_STATE');
        %Clear texts. 
        delete( findobj(gca,'type','text') );
        %Draw selected tracks as green and unselected as red. 
        POS = [];
        c=0;
        hold on
        on_tracks = find( track_matrix(:,t) )';
        for i = on_tracks
            this_track = tracks_curr{i};
            id = t - track_starts(i) + 1;
            POS = this_track(id,3:4);
            %Make clickable text object. 
            text(POS(1),POS(2),num2str(i),'color',c_vec(sel_vec(i)+1,:),'ButtonDownFcn',@textcallback,'UserData',i)
        end
        if CONN_STATE
            group_sel_vec = getappdata(0,'group_sel_vec');
            U = unique( group_sel_vec );
            hlist = findobj(0,'Tag','GroupList');
            list_string = hlist.String;
            if iscell(list_string)
                n_colors = length(hlist.String);
            else
                n_colors = 1;
            end
           
            group_colors = linspecer(n_colors);
            
            div_tracks = find( group_sel_vec > 0 & track_matrix(:,t)' ==1 );
            if( ~isempty(div_tracks))

                for i = div_tracks
                    i
                    this_track = tracks_curr{i};
                    id = t - track_starts(i) + 1;
                    POS = this_track(id,3:4);
                    c = find(U==group_sel_vec(i))-1;
                    text(POS(1),POS(2),num2str(i),'color',group_colors(c,:),'ButtonDownFcn',@textcallback,'UserData',i);
                end
            end            
        end
        
        hold off

        %Add drawspots here so its' called everytime DrawTracks is. 
        if(states.spots)
            DrawSpots
        end
        if(states.cells)
            DrawCells
        end
        
        if PLOT_STATE
            %Update the indicator
            h = getappdata(0,'plot_time_indicator');
            if ~isvalid(h(1))
                return
            end
            plot_signal = getappdata(0,'plot_signal');
            set(h(1),'XData',[t,t]);
            set(h(2),'XData',t);
            t_idx = plot_signal(:,1)==t;
            set(h(2),'YData',plot_signal(t_idx,2));
        end
    end

%Text click call back function. Define so we can use global sel_vec...
    function textcallback(hobj,event)
        %Check status 
        FLAG_STATE = getappdata(0,'FLAG_STATE');
        CONN_STATE = getappdata(0,'CONN_STATE');
        PLOT_STATE = getappdata(0,'PLOT_STATE');
        
        %Get userdata.
        track_id = hobj.UserData
        
        %If we are flagging, then use red to denote flagged tracks. 
        if FLAG_STATE
            
            %Switch status. 
            if( sel_vec(track_id) )
                %Switch off
                sel_vec(track_id) = 0; 
                hobj.Color = c_vec(1,:);
            else
                %Switch on
                sel_vec(track_id) = 1;
                hobj.Color = c_vec(2,:);
            end

            %Final thing to do is output the flagged tracks. 
            flagged.tracks = find(~sel_vec); 
            flagged.spots  = find(~spot_sel_vec);
            flagged.DateString = datestr(clock);
            %Add flagged data to appdata
            setappdata(0,'flagged',flagged);
            
            
        elseif CONN_STATE
            group_sel_vec = getappdata(0,'group_sel_vec');
            
            %Ignore if the current track is flagged. 
            if ~sel_vec(track_id) 
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
                hlist = findobj(0,'Tag','GroupList');   
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
        %Re-drawing cells. 
        
        if(step.spots)
            DrawSpots
        end
        if(step.cells)
            DrawCells
        end
                
        %If plot state is on, plot data associated with this track. 
        if PLOT_STATE

            data_plot = newFigure(8);
            hold on
            %Plot location is wrt MAIN
            %main_pos = MAIN.Position;
            %data_plot.Position =[main_pos(1) main_pos(2)-0.3*screen(4) 0.2*screen(3) 0.2*screen(4)];
        
            
            this_track = tracks_curr{ track_id };
            data =obj.nuc_cyto_data( this_track(:,5));
            sig = cat(1,data.nuc_mean) ./ cat(1,data.cyto_mean);
            title(['Track: ',num2str(track_id)]);
            %Plot the trace. 
            plot(this_track(:,1),sig,'-*','linewidth',2,'color','k');
            
            %Plot a marker / line for the current time. 
            curr_idx = this_track(:,1)==t;

            box off
            ylim([0.5,2])
            yvals=ylim;
            plot_time_indicator(1) = plot([t,t],yvals,'-k');
            plot_time_indicator(2) = plot(t,sig(curr_idx),'*r','markersize',8);
            
            hold off
            setappdata(0,'plot_time_indicator',plot_time_indicator);
            setappdata(0,'plot_signal',[this_track(:,1),sig]);
        end
                
    end

%Spot circle call back function. Upon click...
    function circ_callback(obj,event)
    
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
        
        
        %fit_id
        %spot_sel_vec( fit_id )
    end

%% Drawing spots
    function DrawSpots(object, eventdata)
        MAIN;
        %Remove old spots. 
        delete( findobj(gca,'Type','hggroup') );
        
        if( get(ChBxFilter_hand,'Value'))
            %Get filter value. 
            thresh_val =str2double(get(thresh_val_hand,'string'));
        else
            thresh_val = 0;
        end
        
        %Get cell tracks that are selected. 
        sel_tracks = find(sel_vec);
        
        %Loop over selected tracks. 
        int = [];
        pos = [];
        
        %Looping over cell_track IDs
        for i = sel_tracks
            
            %Look for spot_tracks assigned to this cell. 
            s_tracks_ids = find( spot_track_assignment==i );
            
            if(isempty(s_tracks_ids))
                continue
            end
            
            %Figure out if any tracks exist at time t for cell i. 
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
                        h = viscircles(pos([2,1]),10,'color','g');
                    else
                        h = viscircles(pos([2,1]),10,'color','r');
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
    function DrawCells( object,eventdata)
        
        %Remove contours. 
        delete( findobj(gca,'type','line') ) 
        
        %Pixels of all segemented cells. 
        if ~isfield( DATA.seg(t).frame_obj,'PixelIdxList')
            return
        end
        
        cell_pixels = DATA.seg(t).frame_obj.PixelIdxList;
        MAIN
        %Figure out dimensions. 
        dim = size(DATA.seg(t).frame_obj.centroids{1},2);
        hold on
        if dim==2

            %Convert to [y,x,z]
            for p = 1:length(cell_pixels)
                
                %Get xy location of positive pixels
                [y,x] = ind2sub(DATA.dims(1:2),cell_pixels{p});
                tmp = zeros(DATA.dims(1),DATA.dims(2));
                tmp(cell_pixels{p}) = 1;
                xlims = min(x):max(x);
                ylims = min(y):max(y);
                C = contourc(tmp,[0.5,0.5]);
                [x,y] = C2xyz(C);
                %Find which contour is right. 
                [~,ind] = max( cellfun('length',x) );


                %Plot contour. 
                plot(x{ind},y{ind},'w','linewidth',2);
            end
            
        elseif dim==3
            
            %Project and make convhull? 
            
            disp('missing code')
        end
        
        hold off

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

%% Custom get_planes function. 
function planes = get_planes(reader,Z,C,t)

    for z = 1:Z
        planes(z) = reader.getIndex(z-1,C-1,t-1)+1;
    end


end