

%Default parameters. 

function params = default_params( params )



%List of all default parameters. 
dparams.view_channel=1;
dparams.channel=1;
S  = fieldnames( dparams );

for i = 1:length(S)
    
    %Check if this field exists. 
    if ~isfield(params,S{i})
        params.(S{i}) = dparams.(S{i});
        %Output this default was used. 
        disp(['Using default ',S{i},' with value: ',num2str(dparams.(S{i}))]);
    end
end



end