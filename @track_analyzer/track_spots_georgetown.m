%Tracking segmented objects based on centroids. Requires tracking function from: http://site.physics.georgetown.edu/matlab/code.html
%6-17-18: Updating to make sure it's compatible with 3D tracking. 
%7-12-18: this version is for tracking spot localizations.

function obj = track_spots_georgetown(obj, params, seg_files)


%Default parameters. 
params = default_params(params);

% Image channel
channel_str = ['seg_channel_',pad(num2str(params.seg_channel),2,'left','0')];


%Parse input. Use input file list if provided. 
if(nargin == 2)
    %Get segmentation files
    seg_files = obj.get_frame_files;
end

debug = 0;
%% Collect all position/time data to pass to georgetown code. 

%Loop over seg_files. Also, get all fit data and store into results. 'ids'
%will correspond to the index of the results structure. 
xyzt = [];
fits_all = [];
count = 0;
for i = 1:length(seg_files)
    load(seg_files{i},'frame_obj');
    
    if isfield(frame_obj.(channel_str),params.fit_type)
        these_fits = frame_obj.(channel_str).(params.fit_type);
    else
        these_fits = frame_obj.(params.fit_type);
    end
    
    if(~isempty(these_fits))
        %Position vector. comes as y,x, (z)
        pos = cat(1,these_fits.pos);
        %Switch to x,y (z)
        if(size(pos,2)==3)
            pos = pos(:,[2,1,3]);
            pos_dims = 3;
        elseif(size(pos,2)==2)
            pos = pos(:,[2,1]);
            pos_dims = 2;
        end  
        
        %Filter out dim spots. 
        int = cat(1,these_fits.sum_int);
        snr = cat(1,these_fits.snr);
        
        sel = int >= params.min_intensity & snr >= params.min_snr;
        
        % Sizes could be missing.
        if ~isfield(these_fits,'size')
            warning('No sizes on these fits.');
        else
            sizes = cat(1,these_fits.size);
            sel = sel & sizes >= params.size_range(1) & sizes <= params.size_range(2);
        end
        
        pos = pos(sel,:);
        
        
        %Create ids for all these localizations. 
        ids = [1:size(pos,1)] + count;
        frames = i.*ones(size(pos,1),1);
        %Append to spot matrix ( position, fit_id, frames )
        xyzt = [xyzt; pos, ids', frames ];
        count = size(xyzt,1);
        %Append all_fits (use selection of bright spots)
        fits_all = [fits_all,these_fits(sel)];
    end
end

%Add fits all to obj data. 
obj.results = fits_all;


%Run georgetown tracking (only if multiple time-frames...
if length(seg_files)>1
    
    %set dim parameter. 
    params.dim=pos_dims;
    tracks = track(xyzt,params.max_disp,params);
    
else
    
    tracks = [xyzt, [1:size(xyzt,1)]'];
    
end

if isempty(tracks)
    warning('No tracks found')
    obj.spot_tracks=[];
    return
end


%Reorganize data
ids= unique(tracks(:,end));
spot_tracks = cell( length(ids),1);


%New columns
d=length(tracks(1,:));
new_columns = [d-1,d-2,1:pos_dims];
bad_sum=0;
for i = 1:length(ids)
    sel = tracks(:,end)==ids(i);
    this_spot_track = tracks(sel,new_columns);

    %Find which cell track associates with this spot-track. 
    result_ids = this_spot_track(:,2);
    cell_ids = cat(1,obj.results(result_ids).cell_id);
    frames = this_spot_track(:,1);
    %Find track with matching frame/cell_id combo. 
        v = [frames(1),cell_ids(1)];
        %Loop over cell tracks. Find when track == v. 
        for t = 1:length(obj.tracks)
            cell_frames = obj.tracks{t}(:,1);
            cell_ids    = obj.tracks{t}(:,2);
            ID = find(cell_frames==v(1) & cell_ids==v(2));
            if(~isempty(ID))
                this_track_ID = t;
                break
            end
        end
        
    if(isempty(ID))
        disp('missing track');
        ID = 0;
        bad_sum=bad_sum+1;
    else
        spot_tracks{i} = [this_spot_track(:,1:2),this_track_ID.*ones(size(this_spot_track,1),1),this_spot_track(:,3:end)];
    end            
end

disp(['Total missed:',num2str(bad_sum),' out of ',num2str(length(ids))])
%Add spot_tracks to obj
sel = cellfun(@(x) ~isempty(x), spot_tracks);
obj.spot_tracks = spot_tracks(sel); 

%Write tracking params to exp_info
obj.exp_info.spot_track_params=params;

end



%% Default parameters. 
function step = default_params( step )

%List of all default parameters. 
dstep.min_intensity=0;
dstep.min_snr=0;
dstep.fit_type='fit';
dstep.size_range=[0,Inf];

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
