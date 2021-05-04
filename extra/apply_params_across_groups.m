% Function for applying a set of parameters across an experimental group.
% For now, this function applies steps/params from ALL channels. 

function apply_params_across_groups( csv, exp, groups )

% Load the main experiment. 
info.csv_file=csv;
info.exp_id=exp;
obj=get_exp(info);
obj.exp_info.params.channel_01.cells.px_size=obj.exp_info.pixel_size;

% Get parameters. 
params  =obj.exp_info.params;
steps   =obj.exp_info.steps;

% Loop over group experiments. 
for i = 1:length(groups)
    
    for j = 1:length(groups{i})
        
        info.exp_id=groups{i}(j);
        this_obj=get_exp(info);
        
        this_obj.exp_info.params=params;
        this_obj.exp_info.steps =steps;
        this_obj.save;
        
    end
end

end


