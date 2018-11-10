%% Division detector: uses nuclear area time-trace to detect when divisions happen. 

%Breaks up track into multiple segments if division occurs. Uses window
%parameter and area-change threshold to see if there's a rapid change in
%nuclear area. 

% Important! Assumes only one division per track. 

function split_tracks = division_detector( obj, params )


%Tracks
tracks = obj.tracks;
n_tracks = length(tracks);

%Counter
c=0;
%Output
split_tracks = {};
%Loop over tracks. 
for i = 1:length(tracks)

   %This track. 
   track = obj.tracks{i};
   these_frames = track(:,1);
   
   data = obj.get_track_data( i );
   area_trace = [data.nuc_area].*obj.exp_info.pixel_size(1)^2;
   
   %Track must be long enough for calculations
   if size(track,1) < params.division_window + 1
       continue
   end
   
   %Detector. 
   del_area = diff_arb(area_trace,params.division_window);
   
   %Find when del_area above threshold. 
   sel = del_area <= params.del_area_threshold;
   
   %If multiple potential drops take point of minimum area
   if( sum(sel) >= 1 )
       %Drop frames
       drop_frames = find(sel);
       min_frames  = drop_frames + params.division_window;
       %Areas at the drops
       these_areas = area_trace( min_frames );
       %Take minimum. 
       [~,min_idx] = min( these_areas );
       %Division event frame. Time of minimum area. 
       division_event_frame = min_frames(min_idx);
       %Cut the track off a few frames before drop in nuc-area. Defined by
       %params.cut_frames_before. 
       end_idx = division_event_frame - params.cut_frames_before;
       %Make sure track exists at this time. 
       if end_idx>0
           %Counter
           c=c+1;
           split_tracks{c} = track(1:end_idx,:);
       end

       %Now Check if the track exists a few frames after division. This
       %delay should be a bit longer...?
       start_idx = division_event_frame + params.cut_frames_after;
       
       if start_idx < size(track,1)
          c=c+1; 
          split_tracks{c} = track(start_idx:end,:);
       end
   else
       %No drop, don't chop track. 
       c=c+1;
       split_tracks{c} = track;
   end
   
end
















end