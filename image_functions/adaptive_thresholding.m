%% Thresholdings.
function BW = adaptive_thresholding(I, J, params)


    p_val = params.percentile;

    %Threshold starting point. 
    params.thresh_start = prctile(J(:),p_val);
        if params.thresh_start == 0
            warning('Threshold too low');
            %Calculate the first non-zero percentile and set as
            %starting percentile. 
            ps = prctile(J(:),[p_val:0.2:100]);
            nonzero = ps > 0;
            params.thresh_start = nonzero(1);

        end

    %Incremental change in threshold. 
    %params.increment = params.thresh_start*0.01;
    %First smooth image. 
    I_sm = imgaussfilt(I,params.I_sm_sigma);
    %Perform iterative thresholding. 
    BW = iterative_thresholding(I_sm, J, params );
end