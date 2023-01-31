function varargout = TOOL_BOX(varargin)
% TOOL_BOX MATLAB code for TOOL_BOX.fig
%      TOOL_BOX, by itself, creates a new TOOL_BOX or raises the existing
%      singleton*.
%
%      H = TOOL_BOX returns the handle to a new TOOL_BOX or the handle to
%      the existing singleton*.
%
%      TOOL_BOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOOL_BOX.M with the given input arguments.
%
%      TOOL_BOX('Property','Value',...) creates a new TOOL_BOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TOOL_BOX_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TOOL_BOX_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TOOL_BOX

% Last Modified by GUIDE v2.5 24-Jan-2020 13:12:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TOOL_BOX_OpeningFcn, ...
                   'gui_OutputFcn',  @TOOL_BOX_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before TOOL_BOX is made visible.
function TOOL_BOX_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TOOL_BOX (see VARARGIN)

% Choose default command line output for TOOL_BOX
handles.output = hObject;
handles.draw_cells = @draw_tracks_button_Callback;
handles.step = struct;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TOOL_BOX wait for user response (see UIRESUME)
% uiwait(handles.TOOL_BOX);
%divisions = getappdata(0,'divisions');
%hObject.Position = [605.8571 15 33.1429 34.1875];


% --- Outputs from this function are returned to the command line.
function varargout = TOOL_BOX_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in flag_tracks_button.
function flag_tracks_button_Callback(hObject, eventdata, handles)
% hObject    handle to flag_tracks_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button_state = get(hObject,'Value');

states =getappdata(0,'states');

%Telling what the state of the system is. 
if(button_state)
    %Change background
    hObject.BackgroundColor = [0,1,0];
    
    %Set state on for flag
    states.flag_tracks=1;
       
    %Turn off other states that conflict. 
    states.group=0;
    states.flag_cells=0;
    states.flag_spots=0;
    h(1) = findobj('Tag','CONNECT_BUTTON');
    h(2) = findobj('Tag','flag_cells_button');
    h(3) = handles.flag_spots_button;

    set(h,'BackgroundColor', [1,1,1]);
    set(h,'Value',0);
    
else
    hObject.BackgroundColor = [1,1,1];
    
    %Set state
    states.flag_tracks=0;
end

%Set data
setappdata(0,'states',states);
drawnow

% --- Executes on button press in flag_cells_button.
function flag_cells_button_Callback(hObject, eventdata, handles)
% hObject    handle to flag_cells_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button_state = get(hObject,'Value');
states =getappdata(0,'states');

%Telling what the state of the system is. 
if(button_state)
    %Change background
    hObject.BackgroundColor = [0,1,0];
    
    %Set state on for flag
    states.flag_cells=1;
    %Turn off other states that conflict. 
    states.group=0;
    states.flag_tracks=0;
    states.flag_spots=0;
    h(1) = handles.CONNECT_BUTTON;
    h(2) = handles.flag_tracks_button;
    h(3) = handles.flag_spots_button;
    set(h,'BackgroundColor', [1,1,1]);
    set(h,'Value',0);
    
else
    hObject.BackgroundColor = [1,1,1];
    
    %Set state
    states.flag_cells=0;
end

%Set data
setappdata(0,'states',states);
drawnow;

% --- Executes on button press in flag_spots_button.
function flag_spots_button_Callback(hObject, eventdata, handles)
% hObject    handle to flag_spots_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flag_spots_button
button_state = get(hObject,'Value');
states =getappdata(0,'states');

%Telling what the state of the system is. 
if(button_state)
    %Change background
    hObject.BackgroundColor = [0,1,0];
    
    %Set state on for flag
    states.flag_spots=1;
    %Turn off other states that conflict. 
    states.group=0;
    states.flag_tracks=0;
    states.flag_cells=0;
    h(1) = handles.CONNECT_BUTTON;
    h(2) = handles.flag_tracks_button;
    h(3) = handles.flag_cells_button;
    set(h,'BackgroundColor', [1,1,1]);
    set(h,'Value',0);
    
else
    hObject.BackgroundColor = [1,1,1];
    
    %Set state
    states.flag_spots=0;
end

%Set data
setappdata(0,'states',states);
drawnow;


% --- Executes on button press in CONNECT_BUTTON.
function CONNECT_BUTTON_Callback(hObject, eventdata, handles)
% hObject    handle to CONNECT_BUTTON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button_state = get(hObject,'Value');
states =getappdata(0,'states');

if(button_state)
    %Change background
    hObject.BackgroundColor = [0,1,0];
    %Turn on connector
    states.group=1;
    
    
    states.flag_tracks=0;
    states.flag_cells=0;
    states.flag_spots=0;
    h(1) = handles.flag_tracks_button;
    h(2) = handles.flag_cells_button;
    h(3) = handles.flag_spots_button;

    set(h,'BackgroundColor', [1,1,1]);
    set(h,'Value',0);
else
    hObject.BackgroundColor = [1,1,1];
    
    %Set state
    states.group=0;
end

%Set data
setappdata(0,'states',states);
drawnow;

% --- Executes on button press in SAVE.
function SAVE_Callback(hObject, eventdata, handles)
% hObject    handle to SAVE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Define channel str. 
step = getappdata(0,'step');
channel_str = ['seg_channel_',pad(num2str(step.channel),2,'left','0')];
disp('Saving data....')
%Get track_obj
track_obj = getappdata(0,'track_obj');
%Get flaggs. 
exp_info.flagged = getappdata(0,'flagged');
%Get dividers 
exp_info.groups = getappdata(0,'groups');
track_obj = track_obj.update_exp_info(exp_info);


%Export the drawn cells. 
drawn_cells = getappdata(0,'drawn_cells');

%Img dimensions. 
frame_files = track_obj.get_frame_files();

if(~exist(track_obj.exp_info.nuc_seg_dir,'dir'))
    mkdir(track_obj.exp_info.nuc_seg_dir)
end

%Loop over frames where there was a new drawn cell. 
if ~isempty( drawn_cells )
    disp('Exporting cells...');
    drawn_cells_frames = unique( [drawn_cells.t] );

    for f = drawn_cells_frames

        %Find cells drawn in this frame. 
        idx = find( [drawn_cells.t] == f );

        %Load frame file.
        if exist(frame_files{f}, 'file')
            
            
            
            F = load(frame_files{f});
            [Y,X,Z] = size(F.frame_obj.(channel_str).BW);    

            % Define the dimensionality. 
            dims = length(F.frame_obj.(channel_str).centroids{1});
        
        else
            F=struct('frame_obj',[]);
            F.frame_obj=struct(channel_str,[]);
            Y=track_obj.exp_info.img_size(1);
            X=track_obj.exp_info.img_size(2);
            Z=track_obj.exp_info.z_planes;
            dims=2;
            F.frame_obj.(channel_str)=struct('BW',false([Y,X]),'PixelIdxList',cell(1),'centroids',cell(1),'contours',cell(1));
            frame_files{f} = [track_obj.exp_info.nuc_seg_dir,'frame_',sprintf('%04d',f),'.mat'];
            
        end
        
        %Loop over cells.  
        for c = idx
            % Make a smooth contour around drawn cell. 
            %C = contourc(drawn_cells(c).mask,[0.5,0.5]);
            Cs = drawn_cells(c).vector;
            
            
            % Smoothing coordinates.
            boundaryLength = size(Cs,1);
            
            if boundaryLength < 2
                continue
            end
            
            t = 1 : boundaryLength;
            tq = linspace(1, boundaryLength, 200);
            px = pchip(t, Cs(:,1), tq);
            py = pchip(t, Cs(:,2), tq);
            
            this_mask = poly2mask(px,py,Y,X);
            
            %Check dimensions. 
            [thisY,thisX] = size(this_mask);
            if( thisY~=Y || thisX~=X )
                error('mis-matched dimensions?');
            end
                
            % Process differently depending on 2D/3D.
            if dims == 2
                
                full_mask = this_mask; 
                idx = find(full_mask);
                [y,x] = ind2sub([Y,X],idx);
                ctr = [mean(x),mean(y)];

            elseif dims == 3
            
                full_mask = repmat(this_mask,[1,1,Z]);
                idx = find(full_mask);
                [y,x,z] = ind2sub([Y,X,Z],idx);
                ctr = [mean(x),mean(y),mean(z)];

            end
            
            %Append frame_obj with new cell. 
            F.frame_obj.(channel_str).PixelIdxList{end+1} = find(full_mask);
            F.frame_obj.(channel_str).centroids{end+1} = ctr;
            F.frame_obj.(channel_str).contours{end+1} = [px(:),py(:)];
        end

        % Remove duplicates that were drawn repeatedly. This is needed when
        % you do intermediate saving. 
        all_ctr = round(cell2mat(F.frame_obj.(channel_str).centroids'));                
        D_mat = pdist2(all_ctr,all_ctr);
        sel_mat = D_mat==0 & ~eye(size(all_ctr,1));
        [i,~] = find(tril(sel_mat));
        sel_indices = setdiff([1:size(all_ctr,1)],i);
        F.frame_obj.(channel_str).PixelIdxList =  F.frame_obj.(channel_str).PixelIdxList(sel_indices);
        F.frame_obj.(channel_str).centroids = F.frame_obj.(channel_str).centroids(sel_indices);
        F.frame_obj.(channel_str).contours  = F.frame_obj.(channel_str).contours(sel_indices);
        
        % Update the BW. 
        F.frame_obj.(channel_str).BW = zeros(size(F.frame_obj.(channel_str).BW));
        px = cat(1,F.frame_obj.(channel_str).PixelIdxList{:});
        F.frame_obj.(channel_str).BW(px) = 1;
        
        %Save frame_obj. 
        frame_obj = F.frame_obj;
        save(frame_files{f},'frame_obj');
    end
    exp_info.drawn_cells = drawn_cells;
    exp_info.drawn_cells_frames = unique(drawn_cells_frames);
    track_obj = track_obj.update_exp_info(exp_info);
    try
        params.max_dist = track_obj.exp_info.max_dist;
    catch
        warning('No max dist found. using 1');
        params.max_dist=1;
    end
    
    %Now we re-track to include the new cells. Keep in mind, the gui
    %doesn't update with these new tracks!!
    params.seg_channel=step.channel;
    track_obj = track_obj.track_cells(params);    
end

%Save it. 
track_obj.save;
disp('Finished');
drawnow;

% --- Executes on selection change in GroupList.
function GroupList_Callback(hObject, eventdata, handles)
% hObject    handle to GroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns GroupList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GroupList
curr_entry = hObject.Value;
S = hObject.String;

if iscell(S)
    
    this_string = S(curr_entry);
    split_str = split(this_string,'_');
    curr_id = str2num( split_str{2} );
else
    this_string = S;
    if ~isempty(this_string)
        split_str = split(this_string,'_');
        curr_id = str2num(split_str{2} );
    else
        error('Add track please.');
    end
end
setappdata(0,'group_id',curr_id);
drawnow;

% --- Executes during object creation, after setting all properties.
function GroupList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Get any pre-existing groups in the app data 
D = getappdata(0,'groups');
L = length(D);

if(L==1)
    %Check if empty. 
    if(isempty(D(L).group_id))
        %Set string to empty.
        S = [];
    else
        id = D(L).group_id;
        S = ['Group_',num2str(id)];
    end
elseif(L>1)
    %Make cell_string list. 
    S = cell(L,1);
    for i = 1:L
        id = D(i).group_id;
        S{i} = ['Group_',num2str(id)];
    end
end

hObject.String = S;
curr_entry = hObject.Value;
setappdata(0,'group_id',D(curr_entry).group_id);

% --- Executes on button press in ADD_GROUP.
function ADD_GROUP_Callback(hObject, eventdata, handles)
% hObject    handle to ADD_GROUP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

list_handle = handles.GroupList;

OG_string = list_handle.String;

if(isempty(OG_string))
    N = 0;
    new_string = ['Group_1'];
    max_val=0;
elseif ~iscell(OG_string)
    N=1;
    %Separate og string
    split_str = split(OG_string,'_');
    t_vals = str2num(split_str{end});
    max_val = t_vals;
    new_string{1} = OG_string;
    new_string{2} = ['Group_',num2str(max_val+1)];
else
 
    N = length(OG_string);
    for i = 1:N
        s = OG_string{i};
        split_str = split(s,'_');
        t_vals(i) = str2num(split_str{end});
    end
    max_val = max(t_vals);  
    new_string = OG_string;
    new_string{N+1} = ['Group_',num2str(max_val+1)];
end

%Add another entry list
list_handle.String = new_string;
%Set selection to last string
list_handle.Value = N+1;

setappdata(0,'group_id',max_val + 1);
%Need to add entry into groups structure
G = getappdata(0,'groups');
if isempty(G)
    G = group_structure(1);
end

if isempty(G(1).group_id)
    G(1).group_id=max_val+1;
    G(1).cell_tracks = [];
else
    L = length(G);
    %G(L+1) = group_structure(1);
    G(L+1).cell_tracks = [];
    G(L+1).group_id = max_val + 1;
end
setappdata(0,'groups',G);
drawnow


% --- Executes on button press in REMOVE_GROUP.
function REMOVE_GROUP_Callback(hObject, eventdata, handles)
% hObject    handle to REMOVE_GROUP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%list_handle = findobj('Tag','GroupList');
list_handle = handles.GroupList;
curr_entry = list_handle.Value;

if(~iscell( list_handle.String))
    OG_string = cellstr( list_handle.String);
else
    OG_string =  list_handle.String;
end

%Ignore currently selected string
N = length(OG_string);
if( N==1)
    new_string = [];
    list_handle.String = [];
    list_handle.Value = 1;
    setappdata(0,'group_id',[]);
    setappdata(0,'groups',group_structure(1));
    V = getappdata(0,'group_sel_vec');
    V = zeros(size(V));
    setappdata(0,'group_sel_vec',V);
else
    
    %Locate bad entry. 
    split_str = split( OG_string( curr_entry ),'_' );
    rm_id = str2num( split_str{2} );
    
    sel = [1:N] ~= curr_entry;
    new_string = OG_string(sel);

    selected_string = new_string(N-1);
    split_str = split(selected_string,'_');
    curr_id = str2num(split_str{2});

    %Update object and data. 
    list_handle.String = new_string;
    list_handle.Value = N-1;
    setappdata(0,'group_id',curr_id);
    
    %Get group data. 
    G = getappdata(0,'groups');
    ids = [G.group_id];
    sel = ids ~= rm_id;
    G = G(sel);
    if isempty(G)
        G = group_structure(1);
    end
    setappdata(0,'groups',G);

    %Adjust group_sel_vec
    V = getappdata(0,'group_sel_vec');
    set_these_to_zero = find(V==rm_id);
    V(set_these_to_zero) = 0;
    setappdata(0,'group_sel_vec',V);
    
end
drawnow;


% --- Executes during object creation, after setting all properties.
function CONNECT_BUTTON_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CONNECT_BUTTON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hObject.Value = 0;
states= getappdata(0,'states');
states.group=0;
setappdata(0,'states',states);


% --- Executes during object creation, after setting all properties.
function flag_tracks_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flag_tracks_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hObject.Value = 0;
states= getappdata(0,'states');
states.flag_tracks=0;
setappdata(0,'states',states);



% --- Executes on button press in PLOT_DATA_BUTTON.
function PLOT_DATA_BUTTON_Callback(hObject, eventdata, handles)
% hObject    handle to PLOT_DATA_BUTTON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PLOT_DATA_BUTTON
button_state = get(hObject,'Value');
states=getappdata(0,'states');

if(button_state)
    %Change background
    hObject.BackgroundColor = [0,1,0];
    %Turn on connector
    states.plot_data=1;
else
    hObject.BackgroundColor = [1,1,1];
    
    %Set state
    states.plot_data=0;
end
setappdata(0,'states',states);
drawnow;


% --- Executes during object creation, after setting all properties.
function PLOT_DATA_BUTTON_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PLOT_DATA_BUTTON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hObject.Value = 0;
states=getappdata(0,'states');
states.plot_data=0;
setappdata(0,'states',states);



% --- Executes on button press in cell_division_mark_button.
function cell_division_mark_button_Callback(hObject, eventdata, handles)
% hObject    handle to cell_division_mark_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function just marks the current time point. Useful for difining
% mitotic time line. Fills in marked_frame field of group data. 
G = getappdata(0,'groups');
list_handle = handles.GroupList;
curr_entry = list_handle.Value;

%For the selected group, mark the 'division' time as the current time
%avlue. 
H = findobj('Tag','Main');
t = getappdata(0,'t');
G(curr_entry).marked_frame = t;
disp('Frame marked for group')
[G.marked_frame]
setappdata(0,'groups',G);
drawnow;


% --- Executes on button press in tracks_check_box.
function tracks_check_box_Callback(hObject, eventdata, handles)
% hObject    handle to tracks_check_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tracks_check_box
val=hObject.Value;
s = getappdata(0,'states');
if(val)
    %Change state. 
    s.tracks=1;
else
    s.tracks=0;
end
setappdata(0,'states',s);
drawnow;


% --- Executes on button press in cells_check_box.
function cells_check_box_Callback(hObject, eventdata, handles)
% hObject    handle to cells_check_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tracks_check_box
val=hObject.Value;
s = getappdata(0,'states');
if(val)
    %Change state. 
    s.cells=1;
else
    s.cells=0;
end
setappdata(0,'states',s);
drawnow;



% --- Executes on button press in spots_check_box.
function spots_check_box_Callback(hObject, eventdata, handles)
% hObject    handle to spots_check_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tracks_check_box
val=hObject.Value;
s = getappdata(0,'states');
if(val)
    %Change state. 
    s.spots=1;
else
    s.spots=0;
end
setappdata(0,'states',s);
drawnow;



% --- Executes during object creation, after setting all properties.
function tracks_check_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tracks_check_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
s = getappdata(0,'states');
hObject.Value=s.tracks;


% --- Executes during object creation, after setting all properties.
function cells_check_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cells_check_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
s = getappdata(0,'states');
hObject.Value=s.cells;

% --- Executes during object creation, after setting all properties.
function spots_check_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spots_check_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

s = getappdata(0,'states');
hObject.Value=s.spots;

% --- Executes on scroll wheel click while the figure is in focus.
function TOOL_BOX_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to TOOL_BOX (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
H = findobj('Tag','Main');
figure(H);


% --- Executes on button press in draw_tracks_button.
function draw_tracks_button_Callback(hObject, eventdata, handles)
% hObject    handle to draw_tracks_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of draw_tracks_button

H = findobj('Tag','Main');
figure(H);
% SET axis to the image area?
this_roi = impoly;
disp(['Captured cell.'])
%Get the current data (frame, other drawn cells), 
t = getappdata(0,'t');
drawn_cells = getappdata(0,'drawn_cells');
%if empty, make structure
if isempty(drawn_cells )
    drawn_cells = struct('mask',[],'t',[],'vector',[]);
    n=0;
else
    n = length(drawn_cells);
end    

%Prepare data and append to newly drawn cells. 
new_cell.mask = this_roi.createMask;
new_cell.t = t;
new_cell.vector = this_roi.getPosition;

drawn_cells(n+1) = new_cell;
setappdata(0,'drawn_cells',drawn_cells);
drawnow;


% --- Executes on button press in reset_flags_button.
function reset_flags_button_Callback(hObject, eventdata, handles)
% hObject    handle to reset_flags_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
H=findobj(0,'Tag','Main');
g = guidata(H);
g.reset_flags();
drawnow;

