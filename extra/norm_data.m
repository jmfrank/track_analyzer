%Simple min / max normalization. 

function out = norm_data(in)


min_val = min(in);
max_val = max(in);

out = ( in - min_val ) ./ (max_val - min_val );


end