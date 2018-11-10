
%Generate an empty structure array with the fields used for nascent RNA
%spot tracking. Structure contains all parameters for full 3D guassian fit.
%Fitting functions add 'NaN' to fields if that parameter wasn't used. 
function empty_structure = gen_fit_struct( N )

%Fields
empty_structure(N) = struct('sigma',[],'pos',[],'int',[],...
    'gauss_bg',[],'real_bg',[],'dist_2_periphery',[],'cell_id',[],'resnorm',[],'raw_sum',[],'corrcoef',[]);

end