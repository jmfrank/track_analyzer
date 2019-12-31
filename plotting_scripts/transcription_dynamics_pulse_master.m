%% Master script for looking at nascent transcription dynamics versus time. 

% Treats each group as a separate experiment that happen in a time
% sequence. For example, a pre- and post- treatment of a set of cells. 


function L_group= transcription_dynamics_pulse_master( info, groups, labels, step, params, c_vec )



%Grouping of pulse durations
L_group = cell( length(groups),1);

%Loop over groups
for g = 1:length(groups)
    
    %Experimental ids for group 'g'
    exp_ids = groups{g};
    
    %Figure out how many frames
    info.exp_id=exp_ids(1);
    track_obj=get_exp(info);
    
    time_vec = track_obj.exp_info.time_series*step.time_scale;
    n_frames = length(time_vec);

    max_frame = 0;
    cell_count = 0;
    
    %Vector of pulse durations
    L = [];
    
    %Loop over experiments
    for i = 1:length(exp_ids)
        exp_ids(i);
        %Load experiment. 
        info.exp_id=exp_ids(i);
        track_obj=get_exp(info);
       
        %Track Ids
        indices = find( cellfun('length',track_obj.tracks) >= params.min_length);
                
        %Check for flagged tracks. 
        if(isfield(track_obj.exp_info,'flagged'))
            flagged = track_obj.get_flagged_tracks;
            indices = setdiff(indices,flagged);           
        end
        
        %Create cell array for each time
        %point, that tells us how many cells there are, and if there's at
        %least one nascent spot detected at that time. 
        track_data = cell(n_frames,1);
        
        %First loop over cell tracks. Fill up track_data cells. 
        for j = 1:length(indices)
            this_track = indices(j);
            frames_present = track_obj.tracks{this_track}(:,1);
            cell_ids = track_obj.tracks{this_track}(:,2);
            
            %Loop over frames, add cell_id (identifier) and zero for
            %tracking spots at each time and cell. 
            for t = 1:length(frames_present)
                this_frame = frames_present(t);
                track_data{this_frame} = [track_data{this_frame};[cell_ids(t),0]];
            end
        end
        

        %Now loop over spot_tracks, adding binary answer if cell has at
        %least one spot. 
        for j = 1:length(track_obj.spot_tracks)

            %Length of track
            len = size(track_obj.spot_tracks{j},1);
            
            frames_present = track_obj.spot_tracks{j}(:,1);

            %Get fit ID's
            ids = track_obj.spot_tracks{j}(:,2);

            %The spot_tracking data is stored as: 
            FITS = track_obj.results(ids);

            %Check which cells this track belongs to. 
            cell_id  = cat(1,FITS.cell_id);
            
            
            %Filter out low intensity spots. 
            switch step.sig_type
            
                case 'integrated'
                
                    %Integrated intensity
                    sigma = prod(cat(1,FITS.sigma).^2,2);
                    sig = cat(1,FITS.int).*sigma;
                    
                    if step.FilterIntegrated
                        sel = sig >= params.MinInt;
                    end
                case 'sum_int'
                    sig = cat(1,FITS.sum_int);
                    %Filter out fits below threshold. 
                    if step.FilterRawSum
                        sel = sig >= params.MinSumInt;
                    end
            end
            
            
            sig = sig(sel);
            frames_present = frames_present(sel);
            cell_id = cell_id(sel);
            
            if length(frames_present) > 1
                
                %Check if the intensity threshold breaks up any tracks. 
                t_diff = diff( frames_present )';

                start_tracks = [1,find(t_diff > 1) + 1];
                stop_tracks  = [find(t_diff > 1), length(frames_present)];
                L = [L,stop_tracks - start_tracks + 1];
                
            else
                %Track is only 1 frame long
                L = [L,1];
            end
            
        end
        
        %Count cell tracks
        cell_count = cell_count + length(indices);
    end
    
    L_group{g} = L';
    
end

%Remove outliers
for g = 1:length(L_group)
    sel = L_group{g} <= params.max_duration;
    L_group{g} = L_group{g}(sel);
end


%If output is asked, then don't plot. 
if(nargout==0)
    newFigure(13)

    %violin([1:length(L_group)],L_group,'facecolor',c_vec,'facealpha',0.3)
    plotSpread(L_group,'distributionColors',c_vec)
    boxplot3([1,2],L_group,0.2,zeros(2,3),1)
    hold on

    cellfun(@(x) mean(x), L_group)
    [yes,p] = ttest2( L_group{1},L_group{2});
end


end



