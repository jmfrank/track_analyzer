%% Master script for plotting nuc/cyto ratio over time. 
%Can now handle nuc area trace plotting. 

function [errors,dat_out] = nuc_cyto_plot_master( info, groups, labels,step, params,color_vec )

%Get defaults. 
step = default_step(step);
params = default_params(params);

%Set up subplotting. 
if step.plotNucArea
    SP = 1;
    subplot(2,1,1)
    hold on
    subplot(2,1,2)
    hold on
else
    SP = 0;
end

dat_out = cell(length(groups));

%Loop over groups
for g = 1:length(groups)
    
    %Experimental ids for group 'g'
    exp_ids = groups{g};

    %Empty group variables
    sig_mat = zeros(1,step.FrameRange(end));
    sig_mat_sq = sig_mat;
    counts = sig_mat;
    area_mat = sig_mat;
    area_mat_sq = sig_mat;
    
    
    max_frame = 0;
    cell_count = 0;
    %Loop over experiments
    for i = 1:length(exp_ids)

        info.exp_id = exp_ids(i);
        obj = get_exp(info);
        
        %Track Ids
        indices = find( cellfun(@(x) size(x,1),obj.tracks) >= params.min_length);
        
        %ignore flagged tracks. 
        if isfield( obj.exp_info,'flagged')
           
            bad_ids = obj.exp_info.flagged.tracks;
            indices = setdiff(indices,bad_ids);
            
        end
        
        %Get background. 
        if isempty(obj.background)
            %warning('Missing background data.');
            background_val=0;
        else
            c_idx = find(strcmp(cellstr(strvcat(obj.background.channel)),'YAP'));
            background_val = obj.background(c_idx).val;            
        end

        
        for j = indices(:)'
            
            
            frames_present = obj.tracks{j}(:,1);
            
            if frames_present(1) > params.max_start_frame
                continue
            end
            
            frame_times    = obj.exp_info.time_series(frames_present).*step.TIME_CONVERSION;
                        
            %The spot_tracking data is stored as: 
            data = obj.get_track_data(j);
            

            
            %Ratio
            RATIO = (cat(1,data.nuc_mean)- background_val)./ ( cat(1,data.cyto_mean) - background_val);
            % RATIO = (cat(1,data.nuc_med)-params.mean_background)./ ( cat(1,data.cyto_med) - params.mean_background);
            %Integrated intensity
            sig = RATIO;
            
            %Cut down data
            sel = frames_present <= step.n_frames;
            frames_present = frames_present(sel);
            frame_times = frame_times(sel);
            sig = sig(sel);

                
            %If there are still frames, then continue
            if( isempty( sig ))
                continue
            end

            %Normalize signal to first frame if requested
            if( step.NormSignal )
                sig = sig ./ sig(1);
            end
            
            %Ignore NAN
            if sum( isnan(sig))>0
                continue
            end
            
            %Ignore signals if they pass a max threshold
            if any( sig >= params.max_signal ) | any( sig <= params.min_signal)
                continue
            end
            
            %Ignore cells with low signal. 
            if step.threshold_int
                %Get mean nuc signal. 
                mu_nuc = mean(cat(1,data.nuc_mean));
                if mu_nuc < params.threshold_int
                    continue
                end
            end

            %Area trace.
            if step.plotNucArea
                area_trace = cat(1,data.nuc_area);
                area_trace = area_trace ./ area_trace(1);
                area_trace = area_trace(sel);
                area_mat(frames_present) = area_mat(frames_present) + area_trace';
                area_mat_sq(frames_present) = area_mat_sq(frames_present) + [area_trace.^2]';
                
            end
            
            
            if(step.PlotTraces)
                if(~isfield(step,'alpha'))
                    step.alpha = 0.2;
                end
                if SP
                    subplot(2,1,1)
                end
                plot( frame_times,sig,'linewidth',2,'color',[color_vec(g,:),step.alpha]);
            end
            
            %Output processed data for publishing. 
            dat_out{g}{end+1} = [frames_present, sig];
            
            %Find overlap of frames. 
            [~,idx_signal,idx_sig_mat] = intersect(frames_present,step.FrameRange);
                       
            %Accumulate signal value ('int') over all traces
            sig_mat( idx_sig_mat ) = sig_mat( idx_sig_mat ) + sig(idx_signal)';
            sig_mat_sq( idx_sig_mat ) = sig_mat_sq( idx_sig_mat ) + sig(idx_signal).^2';
            counts( idx_sig_mat )  = counts( idx_sig_mat ) + ones( 1, length(idx_signal));
            cell_count=cell_count+1;
        end
        
        %Determine frames in this experiment
        frame_count = length( obj.exp_info.time_series );
        if( frame_count > max_frame )
            max_frame = frame_count;
        end
        %Count cells
        max_frame
    end
    
    %Append, cell count to labels...
    labels{g} = [labels{g},' n=',num2str(cell_count)];

    %Calculate mean signal for group
    mean_sig = sig_mat ./ counts;
    
    max_frame = min(max_frame, step.FrameRange(end));
    %frames = [0:max_frame-1]*obj.exp_info.frame_time.*step.TIME_CONVERSION;
    frames = obj.exp_info.time_series(1:max_frame)*step.TIME_CONVERSION;
    
    if(step.plot)
        if SP
            subplot(2,1,1)
        end
        main_line(g) = plot(frames(1:max_frame),mean_sig(1:max_frame),'linewidth',5,'color',color_vec(g,:));
        if(step.PlotSTD | step.PlotSEM)
            stdevs = sqrt( ( counts.*sig_mat_sq - sig_mat.^2)./ ( counts.*(counts -1)) );
            if(step.PlotSEM)
                errors = stdevs ./ sqrt( counts );
            else
                errors = stdevs;
            end
            if SP
                subplot(2,1,1);
            end
            errorpatch(frames(1:max_frame),mean_sig(1:max_frame),errors(1:max_frame),color_vec(g,:),0.4);
        end

        
        if SP
            subplot(2,1,2)
            mean_area = area_mat ./ counts;
            plot(frames(1:max_frame),mean_area(1:max_frame),'linewidth',5,'color',color_vec(g,:));
            if(step.PlotSTD | step.PlotSEM)
                stdevs = sqrt( ( counts.*area_mat_sq - area_mat.^2)./ ( counts.*(counts -1)) );
                if(step.PlotSEM)
                    errors = stdevs ./ sqrt( counts );
                else
                    errors = stdevs;
                end
                errorpatch(frames(1:max_frame),mean_area(1:max_frame),errors(1:max_frame),color_vec(g,:),0.4);
            end
        end
    
       
    end

end

if(step.plot)
    if SP
        subplot(2,1,1)
    end
    %Plotting options
    set(gca,'fontsize',20,'linewidth',2)
    main_line
    labels
    l = legend(main_line,labels);
    set(l,'box','off')
    
    if SP
        subplot(2,1,2)
        %Plotting options
        set(gca,'fontsize',20,'linewidth',2);
    end
end



end



%% Check for params / steps. 

%Default parameters. 
function params = default_params( params )
%List of all default parameters. 
dparams.min_signal=0;
dparams.max_signal=20;
dparams.max_start_frame=Inf;

S  = fieldnames( dparams );

for i = 1:length(S)
    
    %Check if this field exists. 
    if ~isfield(params,S{i})
        params.(S{i}) = dparams.(S{i});
        %Output this default was used. 
        disp(['Using default ',S{i},' with value: ',num2str(dparams.(S{i}))]);
    end
end



end

%Default steps. 

function step = default_step( step )

%List of all default steps. 
dstep.plotNucArea=0;
dstep.threshold_int=0;

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
