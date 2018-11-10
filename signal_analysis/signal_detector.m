%% Signal detector function. Sliding window approach to look for rapid
%changes. 
function [sig_vec,N] = signal_detector( frames_present, sig, params )

%Window size
W = params.window_size;
%Length of signal
L = length(sig);
%Threshold of significant change. 
thresh = params.delS_threshold;

%Empty state vector. 
sig_vec = zeros(size(sig));
%Loop over possible windows. 
for i = 1:(L-W)
    
    start_frame = i;
    end_frame   = i+W;
    %Abs magnitude of change over the window. 
    delS = abs( sig(end_frame) - sig(start_frame) );
    %Is this change gr
    if(delS >= thresh)
        sig_vec(start_frame:end_frame) = 1;
    end
end

%Now count up distinct, non-overlapping changes. This doesn't work that
%well for window size of 1, as it would group two consecutive signal
%changes together. 
d = diff(sig_vec);
N = sum( d ~= 0 );



end
    