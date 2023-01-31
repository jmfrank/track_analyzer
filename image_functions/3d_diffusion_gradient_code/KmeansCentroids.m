%%  Written by Bhavna Rajasekaran, PhD
%%  Max Planck Institute for Physics of Complex Systems, Dresden, Germany
%%  Contact: bhavna@pks.mpg.de / bhavnarajasekaran@yahoo.com
%%  Reference: Object segmentation and ground truth in 3D embryonic imaging
%%  Refer to DS_GMM_Kmeans_3dnucleisegmentation.m
%% Plots 3D Centroid Position after Kmeans method on the rendered nuclei (refer to RenderNucleiDS.m)

function KmeansCentroids(C,~)
plot3(C(:,1),C(:,2),C(:,3),'d','markerfacecolor',[0.161  0.337  0.655],'MarkerSize',9);
hold on;
%scatter3(Px(:,1),Px(:,2),Px(:,3),length(Px),IDX,'SizeData',3)
legend('Rendered Nuclei','Centroid DS algorithm(Fused)','GMM Centroid','Kmeans Centroid','Location','NorthEastOutside');
legend('boxoff')
view(-66,24);