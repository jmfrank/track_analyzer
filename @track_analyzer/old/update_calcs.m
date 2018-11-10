%Update n/c ratio and other calcs
function obj = update_calcs(obj, exp_info, params)
    
    %First check if data exists, then calc n/c if not...
    if( obj.nuc_cyto_calc == false )
        obj = obj.nuc_cyto_ratio( exp_info, params);
        track_obj = obj;
        save(exp_info.track_file,'track_obj')
    end
    

end