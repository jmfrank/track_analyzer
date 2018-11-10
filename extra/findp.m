%Find a percentile based on a threshold value. 

function percentile = findp( dist, threshold)


p = [50:0.2:100];
P = prctile(dist,p);

percentile = p(find( P >= threshold, 1));

end