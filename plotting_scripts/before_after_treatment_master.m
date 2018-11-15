%Master script for before/after treatment using nascent transcription cell
%lines. This assumes all experiments are the same image frequency. Specify
%the frames before and after for analyzing. 

function [after, before] = before_after_treatment_master( csv_file, pairs, labels, frames_before, frames_after, dt, step, params, MinSumInt)

%Vector of frame indices before / after. 
frames_before_vec = 1:frames_before;
frames_after_vec  = 1:frames_after;

%Actual times before and after treatment. Use this to make mean fraction
%values. 
times_before = -fliplr(0:dt:(frames_before-1)*dt)';
times_after  = [dt:dt:frames_after*dt]';

sum_before = zeros(size(times_before));
sum_before_sq = zeros(size(times_before));
counts_before = zeros(size(times_before));

sum_after  = zeros(size(times_after));
sum_after_sq  = zeros(size(times_after));
counts_after  = zeros(size(times_after));

%Intensity
INT_before = [];
INT_after  = [];
PULSE_before = [];
PULSE_after  = [];

for i = 1:length(pairs)
    
    %define groups
    groups{1} = pairs{i}(1);
    groups{2} = pairs{i}(2);
    
    params.MinSumInt = MinSumInt(i);
    
    % Send to master script to get fraction of cells. 
    [f_curves,t_vals, INT]  = transcription_dynamics_fraction_master( csv_file, groups, labels, step, params);
    % Get pulse durations. 
    PULSE = transcription_dynamics_pulse_master(csv_file,groups,labels,step,params);
    
    %Get before frames vector. 
    these_frames_before = frames_before-length(f_curves{1})+1:frames_before;
    %Sum up the total vector. 
    [~,idx_data,idx_times_before] = intersect( these_frames_before, frames_before_vec);
    
    %Option to normalize to pre-treatment transcription levels. 
    if step.NormSignal
        %Make transcription relative to mean pre-treatment levels. 
        norm_factor = mean( f_curves{1}(idx_data ) );
    else
        norm_factor = 1;
    end
    data = f_curves{1}(idx_data) ./ norm_factor;
        
    sum_before(idx_times_before) = sum_before(idx_times_before) + data;
    sum_before_sq(idx_times_before) = sum_before_sq(idx_times_before) + data.^2;
    counts_before(idx_times_before) = counts_before(idx_times_before) + ones(length(idx_data),1);
    INT_before = [INT_before; INT{1}];
    PULSE_before = [PULSE_before; PULSE{1}];
    before.curves{i} = data;
    
    
    %Get index wrt 'times_after'
    these_frames_after = 1:length(f_curves{2});
    [~,idx_data,idx_times_after] = intersect( these_frames_after, frames_after_vec);

    %Normalize. 
    data = f_curves{2}(idx_data) ./ norm_factor;

    sum_after(idx_times_after) = sum_after(idx_times_after) + data;
    sum_after_sq(idx_times_after) = sum_after_sq(idx_times_after) + data.^2;
    counts_after(idx_times_after) = counts_after(idx_times_after) + ones(length(idx_data),1);
    INT_after = [INT_after; INT{2}];
    PULSE_after = [PULSE_after; PULSE{2}];
    after.curves{i} = data;
    
end


%Simplify output. 
after.times = times_after;
after.sum = sum_after;
after.sum_sq = sum_after_sq;
after.counts = counts_after;
after.INT  = INT_after;
after.PULSE = PULSE_after;
after.norm = norm_factor;

before.times = times_before;
before.sum = sum_before;
before.sum_sq = sum_before_sq;
before.counts = counts_before;
before.INT  = INT_before;
before.PULSE = PULSE_before;
before.norm = norm_factor;


