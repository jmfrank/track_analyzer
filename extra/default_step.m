

%Default parameters. 

function step = default_step( step )



%List of all default parameters. 
dstep.use_specific_img=0;
dstep.FORCE_ALL_FRAMES=1;
dstep.debug=0;
dstep.subtract_bg=0;
dstep.CLAHE=0;
dstep.MeanFilter=0;
dstep.threshold_by_histogram=0;
dstep.nuc_thresh = 0;
dstep.iterative_thresholding=0;
dstep.channel=1;
dstep.merger=0;

S  = fieldnames( dstep );

for i = 1:length(S)
    
    %Check if this field exists. 
    if ~isfield(step,S{i})
        step.(S{i}) = dstep.(S{i});
        %Output this default was used. 
        disp(['Using default ',S{i},' with value: ',num2str(dstep.(S{i}))]);
    end
end



end