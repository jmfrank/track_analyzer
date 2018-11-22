function  imshow3D_filter( Img, initial_thresh)
%IMSHOW3D displays 3D grayscale or RGB images in slice by slice fashion
%with mouse based slice browsing and window and level adjustment control.
%
% Usage:
% imshow3D ( Image )
% imshow3D ( Image , [] )
% imshow3D ( Image , [LOW HIGH] )
%   
%    Image:      3D image MxNxKxC (K slices of MxN images) C is either 1
%               (for grayscale images) or 3 (for RGB images)  
%    [LOW HIGH]: display range that controls the display intensity range of
%                a grayscale image (default: the widest available range)
%
% Use the scroll bar or mouse scroll wheel to switch between slices. To
% adjust window and level values keep the mouse right button pressed and
% drag the mouse up and down (for level adjustment) or right and left (for
% window adjustment). Window and level adjustment control works only for
% grayscale images.
% 
% "Auto W/L" button adjust the window and level automatically for grayscale
% images
%
% While "Fine Tune" check box is checked the window/level adjustment gets
% 16 times less sensitive to mouse movement, to make it easier to control
% display intensity rang.
%
% Note: The sensitivity of mouse based window and level adjustment is set
% based on the user defined display intensity range; the wider the range
% the more sensitivity to mouse drag.
% 
% Note: IMSHOW3DFULL is a newer version of IMSHOW3D (also available on
% MathWorks) that displays 3D grayscale or RGB images from three
% perpendicular views (i.e. axial, sagittal, and coronal).
% 
%   Example
%   --------
%       % Display an image (MRI example)
%       load mri 
%       Image = squeeze(D); 
%       figure, 
%       imshow3D(Image) 
%
%       % Display the image, adjust the display range
%       figure,
%       imshow3D(Image,[20 100]);
%
%   See also IMSHOW.

%
% - Maysam Shahedi (mshahedi@gmail.com)
% - Released: 1.0.0   Date: 2013/04/15
% - Revision: 1.1.0   Date: 2013/04/19
% - Revision: 1.5.0   Date: 2016/09/22
% 

sno = size(Img,3);  % number of slices
S = round(sno/2);

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

[Rmin Rmax] = WL2R(Win, LevV);


axes('position',[0,0.2,1,0.8]), imshow(squeeze(Img(:,:,S,:)), [Rmin Rmax])

FigPos = get(gcf,'Position');
S_Pos = [50 45 uint16(FigPos(3)-100)+1 20];
Stxt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
Wtxt_Pos = [50 20 60 20];
Wval_Pos = [110 20 60 20];
Ltxt_Pos = [175 20 45 20];
Lval_Pos = [220 20 60 20];
BtnStPnt = uint16(FigPos(3)-250)+1;
if BtnStPnt < 300
    BtnStPnt = 300;
end
Btn_Pos = [BtnStPnt 20 100 20];
ChBx_Pos = [BtnStPnt+110 20 100 20];

%New positions
ChBxFilter_Pos = [BtnStPnt+50 20 100 20];
thresh_text_Pos = [280 20 45 20];
thresh_val_Pos  = [330 20 45 20];


%Estimate threshold step size change. 
vals = Img(:);
increment = 5* (max(vals)-min(vals))/100;
if nargin < 2
    initial_thresh = prctile(vals,90);
end
    



if sno > 1
    shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
else
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
end    
ltxthand = uicontrol('Style', 'text','Position', Ltxt_Pos,'String','Level: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
wtxthand = uicontrol('Style', 'text','Position', Wtxt_Pos,'String','Window: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', WFntSz);
lvalhand = uicontrol('Style', 'edit','Position', Lval_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @WinLevChanged);
wvalhand = uicontrol('Style', 'edit','Position', Wval_Pos,'String',sprintf('%6.0f',Win), 'BackgroundColor', [1 1 1], 'FontSize', WVFntSz,'Callback', @WinLevChanged);
Btnhand = uicontrol('Style', 'pushbutton','Position', Btn_Pos,'String','Auto W/L', 'FontSize', BtnSz, 'Callback' , @AutoAdjust);
ChBxhand = uicontrol('Style', 'checkbox','Position', ChBx_Pos,'String','Fine Tune', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', ChBxSz);
%Adding some boxes to do filtering for peak detection
ChBxFilter_hand = uicontrol('Style','checkbox','Position',ChBxFilter_Pos,'String','Filter','BackgroundColor',[0.8 0.8 0.8], 'FontSize', ChBxSz);
thresh_text_hand = uicontrol('Style', 'text','Position', thresh_text_Pos,'String','Threshold: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
thresh_val_hand = uicontrol('Style', 'edit','Position', thresh_val_Pos,'String',sprintf('%6.0f',initial_thresh), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback',@update_thresh)
ChBxFilter_hand.Value=1;

set (gcf, 'WindowScrollWheelFcn', @mouseScroll);
set (gcf, 'ButtonDownFcn', @mouseClick);
set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
set(gcf,'WindowButtonUpFcn', @mouseRelease)
set(gcf,'ResizeFcn', @figureResized)
set(gcf,'keypressfcn',@fh_kpfcn);
update_thresh();


% -=< Figure resize callback function >=-
    function figureResized(object, eventdata)
        FigPos = get(gcf,'Position');
        S_Pos = [50 45 uint16(FigPos(3)-100)+1 20];
        Stxt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
        BtnStPnt = uint16(FigPos(3)-250)+1;
        if BtnStPnt < 300
            BtnStPnt = 300;
        end
        Btn_Pos = [BtnStPnt 20 100 20];
        ChBx_Pos = [BtnStPnt+110 20 100 20];
        ChBxFilter_Pos = [BtnStPnt 20 100 20];
        
        if sno > 1
            set(shand,'Position', S_Pos);
        end
        set(stxthand,'Position', Stxt_Pos);
        set(ltxthand,'Position', Ltxt_Pos);
        set(wtxthand,'Position', Wtxt_Pos);
        set(lvalhand,'Position', Lval_Pos);
        set(wvalhand,'Position', Wval_Pos);
        set(Btnhand,'Position', Btn_Pos);
        set(ChBxhand,'Position', ChBx_Pos);
        %New stuff
        set(ChBxFilter_hand,'Position',ChBxFilter_Pos);
        set(thresh_text_hand,'Position',thresh_text_Pos);
        set(thresh_val_hand,'Position',thresh_val_Pos);
        
    end

% -=< Slice slider callback function >=-
    function SliceSlider (hObj,event, Img)
        S = round(get(hObj,'Value'));
        
        %Changed to accomodate other plots
        H = findobj(gca,'Type','Image');
        set(H,'cdata',squeeze(Img(:,:,S,:)));
        %set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        
        update_thresh
        
    end

% -=< Mouse scroll wheel callback function >=-
    function mouseScroll (object, eventdata)
        UPDN = eventdata.VerticalScrollCount;
        S = S - UPDN;
        if (S < 1)
            S = 1;
        elseif (S > sno)
            S = sno;
        end
        if sno > 1
            set(shand,'Value',S);
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        %Changed to accomodate other plots
        H = findobj(gca,'Type','Image');
        set(H,'cdata',squeeze(Img(:,:,S,:)))
        
        update_thresh
        %set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
    end

% -=< Mouse button released callback function >=-
    function mouseRelease (object,eventdata)
        set(gcf, 'WindowButtonMotionFcn', '')
    end

% -=< Mouse click callback function >=-
    function mouseClick (object, eventdata)
        MouseStat = get(gcbf, 'SelectionType');
        if (MouseStat(1) == 'a')        %   RIGHT CLICK
            InitialCoord = get(0,'PointerLocation');
            set(gcf, 'WindowButtonMotionFcn', @WinLevAdj);
        end
    end

% -=< Window and level mouse adjustment >=-
    function WinLevAdj(varargin)
        PosDiff = get(0,'PointerLocation') - InitialCoord;

        Win = Win + PosDiff(1) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        LevV = LevV - PosDiff(2) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
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
    
% -=< Arrow key changes thresh value  >=-

    function [] = fh_kpfcn(H,E)          
    % Figure keypressfcn
    
    %Check modifier. 
    if ~isempty(E.Modifier)
        %if shift. 
        if strcmp(E.Modifier{1},'shift')
            inc = 0.1*increment;
        end  
    else
        inc = increment;
    end
    
    switch E.Key
        case 'uparrow'
            

            %Increase threshold. 
            thresh_val =str2double(get(thresh_val_hand,'string'));
            new_val = thresh_val + inc;
            set(thresh_val_hand,'string',num2str(new_val));
            update_thresh();

        case 'downarrow'

            %Decrease threshold. 
            thresh_val =str2double(get(thresh_val_hand,'string'));
            new_val = thresh_val - inc;
            set(thresh_val_hand,'string',num2str(new_val));
            update_thresh();

        otherwise  
    end
    end
    function update_thresh(varargin)
         
        %If the filter checkbox is checked, filter
        to_filter = get(ChBxFilter_hand,'Value');
        
        if(to_filter)
             
            %Clear out old spots
            %H = findobj(gca,'Type','hggroup');
            %delete(H)
 
            %Get threshold
            thresh_val =str2double(get(thresh_val_hand,'string'))
            this_plane = Img(:,:,S);
            p_val = findp(this_plane(:),thresh_val)
             
            %Contours for cell tracking
            H = findobj(gca,'Type','contour');
            delete(H)
             
            hold on
            contour(this_plane,[thresh_val,thresh_val],'color','r','linewidth',2)
            hold off
             
%                 %Get regions
%                 im_threshold = Img >= thresh_val;
%                 stats = regionprops(im_threshold,'Centroid');
% 
%                 centroids = cat(1,stats.Centroid);
%                 centroids = centroids(:,[1,2]);
%                 hold on
%                 H = viscircles(centroids,4*ones(length(stats),1));
%                 hold off          
%             
 
        end
         
    end

    
end