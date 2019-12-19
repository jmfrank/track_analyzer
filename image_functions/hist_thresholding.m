
function BW = hist_thresholding(J, percentile)

    P = prctile(J(:),percentile);
    thrshlevel = P;
    BW  = J >= thrshlevel;

end