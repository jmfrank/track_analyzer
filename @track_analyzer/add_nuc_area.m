%Script for adding nuclear area to nuc/cyto data. Should be faster than
%re-calculating everything. This is for old data where nuc area wasn't
%calculated. 



function obj = add_nuc_area( obj )


%Previously calculating segmentation files%Pre-load all seg_files
seg_files = obj.get_frame_files;


%Loop over frames and add nuc area to data. Should be easy. 
for t = 1:length(seg_files)

    %load data
    load(seg_files{t},'frame_obj');
    %Check if there are any segmented objects
    if ~isfield( frame_obj,'PixelIdxList')
        continue
    end
    
    %Get all areas. 
    A = num2cell( cellfun('length',frame_obj.PixelIdxList) );
    %Set nuc_area field. 
    [frame_obj.data(:).nuc_area] = deal(A{:});
    
    %Now save. 
    save(seg_files{t},'frame_obj');
end



end