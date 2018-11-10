%Create a structure for saving relevant image parameters from LSM files

function obj = save(obj)


%Save myself! Use append in case there are other parts to the track_file. 
track_obj = obj;
%if(exist( obj.exp_info.track_file, 'file'))
%    save(obj.exp_info.track_file,'track_obj','-append');
%else
    save(obj.exp_info.track_file,'track_obj');
%end

end




