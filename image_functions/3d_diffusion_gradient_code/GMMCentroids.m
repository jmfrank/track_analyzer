%%  Written by Bhavna Rajasekaran, PhD
%%  Max Planck Institute for Physics of Complex Systems, Dresden, Germany
%%  Contact: bhavna@pks.mpg.de / bhavnarajasekaran@yahoo.com
%%  Reference: Object segmentation and ground truth in 3D embryonic imaging
%%  Refer to DS_GMM_Kmeans_3dnucleisegmentation.m
%% Plots 3D Centroid Position after GMM method on the rendered nuclei (refer to RenderNucleiDS.m)

function GMMCentroids(model,~,~)
scatter3(model.mu(:,1),model.mu(:,2),model.mu(:,3),'r','filled','SizeData',150);
hold on;
%scatter3(new_cluster_cat(:,1),new_cluster_cat(:,2),new_cluster_cat(:,3),length(new_cluster_cat),idxgm,'SizeData',15)