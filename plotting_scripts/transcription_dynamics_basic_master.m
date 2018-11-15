%% Master script for looking at nascent transcription dynamics versus time. 

% Treats each group as a separate experiment that happen in a time
% sequence. For example, a pre- and post- treatment of a set of cells. 


function transcription_dynamics_basic_master( csv_file, groups, labels, step, params, c_vec )

% figure(1)
% clf
% set(gcf,'color','w')
% subplot(2,1,1)
% std_figure_settings
% subplot(2,1,2);
% std_figure_settings

delT = 0;

%Group of intensities
INT_group = cell(length(groups),1);

%Loop over groups
for g = 1:length(groups)
    
    %Experimental ids for group 'g'
    exp_ids = groups{g};
    
    %Figure out how many frames
    exp_info = get_exp_info_lsm_series(csv_file, exp_ids(1));
    load(exp_info.track_file);
    time_vec = track_obj.exp_info.time_series*step.time_scale + delT;
    n_frames = length(time_vec);
    %Empty group variables
    sig_mat = zeros(1,n_frames);
    sig_mat_sq = sig_mat;
    counts = sig_mat;
    total_counts = sig_mat;
    
    max_frame = 0;
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
        indices = find( cellfun('length',track_obj.tracks) >= params.min_length);
                
        %Check for flagged tracks. 
        if(isfield(track_obj.exp_info,'flagged'))
            
            if ~isempty( track_obj.exp_info.flagged )

                flagged_tracks = track_obj.exp_info.flagged.tracks;

                %Remove flagged tracks from indices. 
                indices = setdiff(indices,flagged_tracks);

                %Now we need to get list of cell_ids / frames so we can ignore
                %spots here. 
                bad_dat = [];
                for t = 1:length(flagged_tracks)
                   bad_dat = [bad_dat; track_obj.tracks{ flagged_tracks(t) }(:,1:2)];
                end
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
            
            %Now filter out any data from flagged tracks. 
            
            
            %Loop over each frame in the track and update 'track_data'
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
%             
%             if(step.PlotTraces)
%                 subplot(2,1,2)
%                 hold on
%                 plot(time_vec(frames_present),sig,'linewidth',1,'color',[c_vec(g,:),0.5]);
%                 hold off
%             end

            %Accumulate signal value ('int') over all traces
            sig_mat( frames_present ) = sig_mat( frames_present ) + sig';
            sig_mat_sq( frames_present ) = sig_mat_sq( frames_present ) + sig.^2';
            %Count the number of 'on' cells. Used for mean intensity, and
            %fraction of transcribing cells. 
            counts( frames_present )  = counts( frames_present ) + 1; 
            
            %Keep track of all localization intensities. 
            INT = [INT, sig'];
        end
        
        %Count cell tracks
        cell_count = cell_count + length(indices);
    end
    
    %Append, cell count to labels...
    labels{g} = [labels{g},' n=',num2str(cell_count)];

    %% New way of calculating fraction. 
        total_counts = cellfun(@(x) size(x,1), track_data);
        on_counts = cellfun(@(x) sum(x(:,2)),track_data);
        %Fraction transcribing
        fraction = on_counts ./ total_counts;
    %subplot(2,1,1)
    hold on
    plot(time_vec,fraction,'-','linewidth',3,'color',[c_vec(g,:),0.6]);
    
    %Smooth or mean fraction. 
    switch step.mean_type
    
        case 'mean'
        
            mean_fraction = mean(fraction);
            plot([time_vec(1),time_vec(end)],[mean_fraction,mean_fraction],'--','linewidth',4,'color',[c_vec(g,:),0.8]);
            
        case 'smooth'
            
            smoo = smooth(fraction);
            plot(time_vec,smoo,'-','linewidth',5,'color',[c_vec(g,:),1]);
        
    end
    hold off
    %Plotting mean signal
    mean_sig = sig_mat ./ counts;
    %subplot(2,1,2)
    %hold on
    %plot(time_vec,mean_sig,'linewidth',5,'color',c_vec(g,:));
    %hold off
    
    %Plotting error bars
%     if(step.PlotSTD | step.PlotSEM)
%         stdevs = sqrt( ( counts.*sig_mat_sq - sig_mat.^2)./ ( counts.*(counts -1)) );
%         if(step.PlotSEM)
%             errors = stdevs ./ sqrt( counts );
%         else
%             errors = stdevs;
%         end
%         %Make patch for errors
%         x = time_vec;
%         uE = mean_sig + errors; lE = mean_sig - errors;
%         yP = [lE'; fliplr(uE)'];
%         xP = [ x' ; fliplr(x)' ] ;
%         try
%             patch( xP, yP, 1, 'FaceColor',c_vec(g,:),'EdgeColor','none',...
%                 'FaceAlpha',0.4);
%         catch
%             disp('patch error')
%         end
       
%    end
    delT = time_vec(end);
    INT_group{g} = INT;
end

%subplot(2,1,1)
%ylim([0,1])

%subplot(2,1,2)
%plotSpread(INT_group,'distributionColors',c_vec)
%hold on
%boxplot2([1:length(INT_group)],INT_group,0.1,c_vec,2)
%[yes,p] = ttest2( INT_group{1},INT_group{2})


end



