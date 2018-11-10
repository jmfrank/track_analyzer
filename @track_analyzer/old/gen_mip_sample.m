%Generate an MSD object with the trackanalyzer object? 


function obj = gen_mip_sample(obj)



    tracks = cellfun(@(x) x(:,[1,3,4]), obj.tracks,'UniformOutput',0);

    msd_obj = msdanalyzer(2,'','');
    msd_obj = msd_obj.addAll(tracks);

    obj.msd_obj = msd_obj;
    
end