%% Filter spot fits. 
function filtered = filter_spot_fits(fits, params)

    sel = [fits.real_bg] < params.max_bg_int;
    fits=fits(sel);
    
    sel = [fits.snr] >= params.min_snr;
    fits=fits(sel);
    
    sel = [fits.size] <= params.max_vol & [fits.size] >= params.min_vol;
    fits=fits(sel);
    
    sum_ints = [fits.sum_int];
    sel = sum_ints >= params.min_sum_int & sum_ints <= params.max_sum_int;
    fits=fits(sel);
    
    sig2bg = ([fits.peak]+[fits.real_bg])./[fits.real_bg];
    sel = sig2bg >= params.min_sig2bg;
    filtered=fits(sel);

end