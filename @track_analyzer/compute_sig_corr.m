function obj = compute_sig_corr(obj, varargin)

%% ADAPTED FROM msdanalyzer velocity correlation. 
%%COMPUTEVCORR Compute velocity autocorrelation.
%
% obj = obj.computeVCorr computes the velocity autocorrelation for all
% the particles trajectories stored in this object. Velocity
% autocorrelation is defined as vc(t) = < v(i+t) x v(i) >, the mean
% being taken over all possible pairs inside a trajectories.
%
% Results are stored in the 'vcorr' field of the returned object.
% The velocity autocorrelation is stored for each particles in a
% cell array, one cell per particle. The array is a double array of
% size N x 4, and is arranged as follow: [dt mean std N ; ...]
% where dt is the delay for the autocorrelation, mean is the mean
% autocorrelation value for this delay, std the standard deviation
% and N the number of points in the average.
%
% obj = obj.computeVCorr(indices) computes the velocity
% autocorrelation only for the particles with the specified
% indices. Use an empty array to take all particles.
%% 


%Input parsing. 
p = inputParser;
p.addRequired('params',@isstruct);
p.addParameter('tracks',obj.tracks,@iscell);
p.addParameter('indices','all');
p.parse(varargin{:})

tracks = p.Results.tracks;
params = p.Results.params;
%Indices depends on what type of input. 
if strcmp(p.Results.indices,'all')
    indices = 1:length(tracks);
else
    indices = p.Results.indices;
end

obj.sigcorr = cell(numel(tracks), 1);


% Get some information before computing...

%Max time delay observed 
max_t = min(params.max_length,max( cellfun(@(x) x(end,1), tracks)));
delays = [1:max_t]';
n_delays = numel(delays);
n_tracks = numel(indices);

%Looping over tracks
fprintf('Computing signal autocorrelation of %d tracks... ', n_tracks);
fprintf('%4d/%4d', 0, n_tracks);
for i = 1 : n_tracks
    fprintf('\b\b\b\b\b\b\b\b\b%4d/%4d', i, n_tracks);
    
    % Holder for mean, std calculations
    sum_vcorr     = zeros(n_delays-1, 1);
    sum_vcorr2    = zeros(n_delays-1, 1);
    n_vcorr       = zeros(n_delays-1, 1);
    
    this_track = indices(i);
    % Get result IDS. 
    ids = tracks{this_track}(:,5);
    %Data
    if isfield( obj.nuc_cyto_data,'nuc_mean')

        dat = obj.nuc_cyto_data( ids );

    else
        if isfield( obj.nuc_cyto_data, 'data')
            channel='data';
        else
            channel='channel_01';
        end
        dat = obj.nuc_cyto_data.(channel)( ids );
    end
    
    %Signal type. 
    switch params.sig_type
        
        case 'mean'
            
            sig = ( cat(1,dat.nuc_mean) - params.mean_background) ./ ( cat(1,dat.cyto_mean) - params.mean_background);
            
        case 'median'
            
            sig = ( cat(1,dat.nuc_med ) - params.mean_background) ./ ( cat(1, dat.cyto_med) - params.mean_background);
            
    end
    
    
    %Check for errors. 
    if(sum( isnan(sig) > 0 ))
        this_track
        pause
    end
    
   
    %Cutoff signal. 
    max_val = min(params.max_length, length(sig));
    
    %Decide which type of correlation. 
    switch params.corr_type
        
        case 'std'
            
            t = tracks{this_track}(:, 1);
            sig = sig(1:max_val);
            t = t(1:max_val);
            V = sig;
            n_detections = size(V, 1);

            % First compute correleation at dt = 0 over all tracks
            Vc0 = mean( sum( V.^2, 2) );

            % Other dts
            for j = 1 : n_detections - 1

                % Delay in frames
                dt = t(j+1:end) - t(j);

                % Determine target delay index in bulk
                [~, index_in_all_delays, ~] = intersect(delays, dt);

                % Velocity correlation in bulk
                lvcorr = sum( repmat(V(j, :), [ (n_detections-j) 1]) .* V(j+1:end, :), 2 );

                % Normalize
                lvcorr = lvcorr ./ Vc0;

                % Store for mean computation
                sum_vcorr(index_in_all_delays)   = sum_vcorr(index_in_all_delays) + lvcorr;
                sum_vcorr2(index_in_all_delays)  = sum_vcorr2(index_in_all_delays) + (lvcorr .^ 2);
                n_vcorr(index_in_all_delays)     = n_vcorr(index_in_all_delays) + 1;
                
            end

            mean_vcorr = sum_vcorr ./ n_vcorr;
            std_vcorr = sqrt( (sum_vcorr2 ./ n_vcorr) - (mean_vcorr .* mean_vcorr)) ;
            vcorrelation = [ delays(1:end-1) mean_vcorr std_vcorr n_vcorr ];
            vcorrelation(1,:) = [0 1 0 n_detections];

            % Store in object field
            obj.sigcorr{i} = vcorrelation;
            
        case 'xcorr'
            max_idx = min(length(sig),33);
            [vcorrelation, lags] = xcorr(sig(1:max_idx),'coeff');
            obj.sigcorr{i} = [lags',vcorrelation];
            
        case 'dev'
            
            %Looks at magnitude of signal change at different times.
            t = tracks{this_track}(:, 1);
            sig = sig(1:max_val);
            t = t(1:max_val);
            V = sig;
            n_detections = size(V, 1);

            % First compute correleation at dt = 0 over all tracks
            Vc0 = mean(  V.^2 );

            % Other dts
            for j = 1 : n_detections - 1

                % Delay in frames
                dt = t(j+1:end) - t(j);

                % Determine target delay index in bulk
                [~, index_in_all_delays, ~] = intersect(delays, dt);

                % Velocity correlation in bulk
                dev =  abs( repmat(V(j, :), [ (n_detections-j) 1]) - V(j+1:end, :) );

                % Don't need to normalize. The theoretical deviation at delT = 0 should be zero.              
                %dev = dev ./ Vc0;

                % Store for mean computation
                sum_vcorr(index_in_all_delays)   = sum_vcorr(index_in_all_delays) + dev;
                sum_vcorr2(index_in_all_delays)  = sum_vcorr2(index_in_all_delays) + (dev .^ 2);
                n_vcorr(index_in_all_delays)     = n_vcorr(index_in_all_delays) + 1;
            
            end

            mean_vcorr = sum_vcorr ./ n_vcorr;
            std_vcorr = sqrt( (sum_vcorr2 ./ n_vcorr) - (mean_vcorr .* mean_vcorr)) ;
            vcorrelation = [ delays(1:end-1) mean_vcorr std_vcorr n_vcorr ];
            %vcorrelation(1,:) = [1 0 0 n_detections];
            %Add in a value at t=0. 
            vcorrelation = [0,0,0,0; vcorrelation];
            
            % Store in object field
            obj.sigcorr{this_track} = vcorrelation;
            
    end
       
    
end
fprintf('\b\b\b\b\b\b\b\b\bDone.\n')

end
