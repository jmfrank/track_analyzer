%% Master script for looking at fraction of transcribing cells over time. Uses current axes to plot.

function [f_curves,t_vals,INT_group] = transcription_dynamics_fraction_master( csv_file, groups, labels, step, params )

%Group of intensities
INT_group = cell(length(groups),1);

%Loop over groups
for g = 1:length(groups)
    
    %Experimental ids for group 'g'
    exp_ids = groups{g};
    
    %Figure out how many frames
    exp_info = get_exp_info_lsm_series(csv_file, exp_ids(1));
    load(exp_info.track_file);
    time_vec = track_obj.exp_info.time_series*step.time_scale;
    %Adjust time-vector so group 1 ends at t=0. 
    if(g==1)
        time_vec = time_vec - time_vec(end);
    elseif(g==2)
        time_vec = time_vec + time_vec(2);
    end
    n_frames = length(time_vec);
    
    
    %Empty group variables
    sig_mat = zeros(1,n_frames);
    sig_mat_sq = sig_mat;
    counts = sig_mat;
    cell_count = 0;
    
    %Intensity vector
    INT = [];
    %Loop over experiments
    for i = 1:length(exp_ids)
        exp_ids(i);
        %Get exp info for cell tracking object. Different csv_file type for lsm
        %series
        exp_info = get_exp_info_lsm_series(csv_file, exp_ids(i));

        %Load track file
        load( exp_info.track_file);
        
        %Track Ids
        indices = find( cellfun(@(x) size(x,1),track_obj.tracks) >= params.min_length);
                
        %Check for flagged tracks. 
        if(isfield(track_obj.exp_info,'flagged'))
            %Might be empty. 
            if ~isempty( track_obj.exp_info.flagged )
                flagged_tracks = track_obj.exp_info.flagged.tracks;
                %Remove flagged tracks from indices. 
                indices = setdiff(indices,flagged_tracks);
            end
            

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
            
            %Loop over each frame in the track and update 'track_data'.
            %Remove spots associated with a flaggged track. 
            for t = 1:length(frames_present)
                this_frame = frames_present(t);
                this_cell  = cell_id(t);
                
                %Getting data for this frame
                dat = track_data{ this_frame };
                %Find which row 'this_cell' is at.
                row = find(dat(:,1)==this_cell);
                
                %For the matched rows, change the second column of dat to
                %1 (this means this cell is on). 
                dat(row,2) = 1;
                
                %Now replace track_data with data at this frame
                track_data{ this_frame } = dat;
                
            end


            %Accumulate signal value ('int') over all traces
            sig_mat( frames_present ) = sig_mat( frames_present ) + sig';
            sig_mat_sq( frames_present ) = sig_mat_sq( frames_present ) + sig.^2';
            %Count the number of 'on' cells. Used for mean intensity, and
            %fraction of transcribing cells. 
            counts( frames_present )  = counts( frames_present ) + 1; 
            
            %Keep track of all localization intensities. 
            INT = [INT; [time_vec(frames_present)',sig]];
        end
        
        %Count cell tracks
        cell_count = cell_count + length(indices);
    end
    
    %Append, cell count to labels...
    labels{g} = [labels{g},' n=',num2str(cell_count)];
    
    %Remove empty frames. 
    sel = ~cellfun('isempty',track_data);
    track_data = track_data(sel);
    time_vec = time_vec(sel);
    
    %Caclulate total number of cells / frame and fraction of 'on' cells. 
    total_counts = cellfun(@(x) size(x,1), track_data);
    on_counts = cellfun(@(x) sum(x(:,2)),track_data);
    %Fraction transcribing
    fraction = on_counts ./ total_counts;
    f_curves{g} = fraction;
    t_vals{g} = time_vec;
    %Intensities of spots
    INT_group{g} = INT;
    
end


end



