%Simple min / max normalization. 

function out = norm_data(in,interval)

if nargin < 2
    interval = [0,1];
end

min_val = min(in(:));
max_val = max(in(:));

out =(interval(2)-interval(1)).*( in - min_val ) ./ (max_val - min_val ) + interval(1);


end