%% Master script for looking at pulse duration of nascent transcription spots. 

%Running this requires using a particle tracking algorithm on the nascent
%spot data. 

function transcription_pulse_duration_master( csv_file, groups, labels, step, params, c_vec )
delT = 0;
%Loop over groups
for g = 1:length(groups)
    
    %Experimental ids for group 'g'
    exp_ids = groups{g};
   
    %%%%%Empty group variables
    
    %Vector of transcription duration (frames?)
    dur_vec = [];
    cell_count = 0;
    
    %Loop over experiments
    for i = 1:length(exp_ids)
        
        exp_ids(i);
        
        %Get exp info for cell tracking object. Different csv_file type for lsm
        %series
        exp_info = get_exp_info_lsm_series(csv_file, exp_ids(i));

        %Load track file
        load( exp_info.track_file);
        
        %Track Ids
        indices = find( cellfun('length',track_obj.spot_tracks) >= params.min_length);
        
        %Create a particle tracking object from fitting data. 
            %All data
            data = track_obj.results;
        
        
        
        %Loop over selected tracks
        for j = indices 
            %Length of track
            len = size(track_obj.spot_tracks{j},1);
            %Frames this track is present
            frames_present = track_obj.spot_tracks{j}(:,1);

            
            %Divide spot_tracks into consecutive 
            %Get fit ID's
            ids = track_obj.spot_tracks{j}(:,2);

            %The spot_tracking data is stored as: 
            FITS = track_obj.results(ids);

            %Get ids for fits with good fit (low resnorm value)
            %fit_ids = cat(1,FITS.resnorm) <= params.resnorm_cutoff;
            %FITS = FITS(fit_ids);
            switch step.sig_type
            
                case 'integrated'
                
                    %Integrated intensity
                    sigma = prod(cat(1,FITS.sigma).^2,2);
                    sig = cat(1,FITS.int).*sigma;
                    
                    if step.FilterIntegrated
                        sel = sig >= params.MinInt;
                        sig = sig(sel);
                        frames_present = frames_present(sel);
                    end
                case 'sum_int'
                    sig = cat(1,FITS.sum_int);
                    %Filter out fits below threshold. 
                    if step.FilterRawSum
                        sel = sig >= params.MinRawSum;
                        sig = sig(sel);
                        frames_present = frames_present(sel);
                    end
                    
                case 'raw_sum'
                    sig = cat(1,FITS.raw_sum);
                    %Filter out fits below threshold. 
                    if step.FilterRawSum
                        sel = sig >= params.MinRawSum;
                        sig = sig(sel);
                        frames_present = frames_present(sel);
                    end
            end
            
            %Get rid of empty tracks. 
            if(isempty(FITS))
                continue
            end
            
            %Normalize signal to first frame if requested
            if( step.NormSignal )
                sig = sig ./ sig(1);
            end

            %Count the number of 'on' cells. Used for mean intensity, and
            %fraction of transcribing cells. 
            counts( frames_present )  = counts( frames_present ) + 1; 
            %Total number of cells at this time. 
            all_frames_present = track_obj.tracks{j}(:,1);
            total_counts( all_frames_present ) = total_counts( all_frames_present ) + 1;
        end
        
        %Count cells
        cell_count = cell_count + length(indices);
    end
    
    %Append, cell count to labels...
    labels{g} = [labels{g},' n=',num2str(cell_count)];

    %Fraction transcribing
    fraction = counts ./ total_counts .* 100;
    hold on
    h = plot(time_vec,fraction,'linewidth',5,'color',c_vec(g,:));
    sel = ~isnan(fraction);
    if(step.mean_val)
        mean_fraction = mean(fraction(sel));
        plot([time_vec(1),time_vec(end)],[mean_fraction,mean_fraction],'--','linewidth',4,'color','k'); %[c_vec(g,:),0.8]);
    end
    hold off

    f_curves{g} = fraction;
    t_vals{g} = time_vec;
end

ylim([0,1])


end


