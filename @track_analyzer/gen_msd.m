%Generate an MSD object with the trackanalyzer object? 


function msd_obj = gen_msd(obj)



    tracks = cellfun(@(x) x(:,[1,3,4]), obj.tracks,'UniformOutput',0);

    msd_obj = msdanalyzer(2,'','');
    msd_obj = msd_obj.addAll(tracks);
    
end