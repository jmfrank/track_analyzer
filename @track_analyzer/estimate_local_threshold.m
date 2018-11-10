

%Testing blockproc. 

function T = estimate_local_threshold( block_struct, hist_tresh, thrshlevel )

    p_vals = [40:0.5:100];
    P = prctile(block_struct.data,p_vals);
    dP = diff(P);
    
    [pks,locs,w,p] = findpeaks(dP);
    
    %Peaks above prominance 
    peak_idx = find(p > params.d_hist_threshold);
    if( ~isempty(peak_idx) )
       %Obvious peak found. Use this one. 
       ptile_idx = locs(peak_idx(1));
       thrshlevel = P(ptile_idx);

    elseif( isempty(peak_idx) )

       %All else fails, try finding intersection
       k = find(dP>=hist_tresh);
       if( ~isempty(k) )

            thrshlevel = P(k(1));
       end

   end

   %Increase 10%
   thrshlevel = 1.1*thrshlevel;

end