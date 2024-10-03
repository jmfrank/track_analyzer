%% FIXING old style data structure. 


load(F{1});

p_list = frame_obj.seg_channel_02.cells.PixelIdxList;

%%

img_size = obj.exp_info.img_size;

BW = false([img_size,obj.exp_info.z_planes]);
BW( cat(1,p_list{:}) ) = 1;

%% calculate stats object. 

stats = table2struct(regionprops3(BW, 'VoxelIdxList','Centroid','SurfaceArea','EigenVectors','EigenValues','Volume'));
[stats.PixelIdxList] = stats.VoxelIdxList;
stats = rmfield(stats,'VoxelIdxList');

%% Save stats to frame_obj. 

%frame_obj.seg_channel_02.cells = rmfield(frame_obj.seg_channel_02.cells,'PixelIdxList');
frame_obj.seg_channel_02.cells.stats = stats;

save(F{1}, 'frame_obj','-v7.3');

