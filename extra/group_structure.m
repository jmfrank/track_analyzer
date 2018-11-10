%Structure for defining groups of cells. 
function D = group_structure(N)
D(N) = struct('group_id',[],'cell_tracks',[],'marked_frame',[]);
end