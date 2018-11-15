%% Division detector: uses nuclear area time-trace to detect when divisions happen. 

% Breaks up track into multiple segments if division occurs. Uses window
% parameter and area-change threshold to see if there's a rapid change in
% nuclear area.
% 
% Should be able to handle multiple divisions per track. Uses peakfinding
% to look for rapid negative chances in nuclear area

function [split_tracks, src_tracks] = division_detector( obj, params, step )

if nargin < 3
    step = struct;
    step = default_step(step);
end

%Tracks
tracks = obj.tracks;
n_tracks = length(tracks);

%Counter
c=0;
%Output
split_tracks = {};
src_tracks=[];
%Loop over tracks. 
for i = 1:length(tracks)
    
    
    %This track. 
    track = obj.tracks{i};
    these_frames = track(:,1);

    data = obj.get_track_data( i );
    area_trace = [data.nuc_area].*obj.exp_info.pixel_size(1)^2;

    %Track must be long enough for calculations
    if size(track,1) < params.division_window + 1 | size(track,1) < params.min_length
       continue
    end

    %Detector. 
    del_area = diff_arb(smooth(area_trace,3),params.division_window);

    %Find where del-area below threshold. 
    sel = del_area < params.del_area_threshold;

    %Find regions above threshold. 
    [vals,lengths,bi] = RunLength(sel);

    %Number of potential divisions/big area changes. 
    loc = find(vals==1);
    N = sum(vals==1);
   
    %Loop over vals. Break up tracks. 
    start_frame=1;
    
    if N>0

        for n = 1:N

            %Check startframe
            if start_frame >= size(track,1)
                break
            end

            idx = loc(n);

            %Get frame vals for this region. 
            frame_vals = [bi(idx):bi(idx)+lengths(idx)-1]+params.division_window;

            %Find frame of minimum area.
            [~,min_idx] = min( area_trace(frame_vals) );

            %Define the division event frame (minimum area)
            division_event_frame = frame_vals(min_idx);

            %Estimate segment end frame. 
            end_frame = division_event_frame - params.cut_frames_before;

            %Check if that end frame after start_frame. 
            if end_frame > start_frame
               %Counter
               c=c+1;
               split_tracks{c} = track(start_frame:end_frame,:);
               src_tracks(c) = i;
               %Move start_frame. 
               start_frame = division_event_frame + params.cut_frames_after;
               
            else
                %Skip this one. 
                start_frame=division_event_frame+params.cut_frames_after;
                continue
            end      
        end

        %Now clean up last segment. 
        if start_frame < size(track,1)
            end_frame = size(track,1);
            c=c+1;
            split_tracks{c} = track(start_frame:end_frame,:);
            src_tracks(c) = i;
        end
        
        if step.debug
           newFigure(2)
           plot(area_trace,'linewidth',2)
        end
            
    else
       %No drop, don't chop track. 
       c=c+1;
       split_tracks{c} = track;        
       src_tracks(c) = i;
    end
   
end
















end                      