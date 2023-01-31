%%  Written by Bhavna Rajasekaran, PhD
%%  Max Planck Institute for Physics of Complex Systems, Dresden, Germany
%%  Contact: bhavna@pks.mpg.de / bhavnarajasekaran@yahoo.com
%%  Reference: Object segmentation and ground truth in 3D embryonic imaging
%%  Refer to DS_GMM_Kmeans_3dnucleisegmentation.m

%%  Renders 3D nuclei for object under-going post-processing & Centroid postion by DS Algorithm

function RenderNucleiDS(G,x,y,z,freqx,freqy,freqz,X_Centroid_px,Y_Centroid_px,Z_Centroid_px)
%%
figure;
border=2;
[nx,ny,nz]=size(G);

Pxls_yrange=max([1 (freqx(2,1)-border)]):min([(freqx(end,1)+border) ny]);
Pxls_xrange=max([1 (freqy(2,1)-border)]):min([(freqy(end,1)+border) nx]);
Pxls_zrange=max([1 (freqz(2,1)-border)]):min([(freqz(end,1)+border) nz]);

Irange=G(Pxls_xrange,Pxls_yrange,Pxls_zrange);
[X Y Z]  = meshgrid(y(Pxls_yrange),x(Pxls_xrange),z(Pxls_zrange));
Thresh=max(double(Irange(:)))*0.2;
hpatch = patch(isosurface(X,Y,Z,Irange,Thresh));
isonormals(X,Y,Z,Irange,hpatch);
set(hpatch,'FaceColor',[1 0.843137 0],'EdgeColor','none','DiffuseStrength',0.4,'SpecularStrength',0.4,'FaceAlpha',0.3);
set(gca,'FontSize',7);
set(gca,'TickLength',[0 0],'XMinorTick','off','YMinorTick','off','ZMinorTick','off');
set(gca,'LineWidth',1);
xlabel('X ({\mu}m) ','FontSize',14);
ylabel('Y ({\mu}m) ','FontSize',14);
zlabel('Z ({\mu}m)','FontSize',14);
view(3);
daspect([1 1 1]);
camlight left;
%set(gcf,'Renderer','zbuffer')
% axis tight;
set(gcf,'Renderer','opengl');
set(gca,'TickLength',[0 0]);
set(gca,'FontSize',18);
set(gca,'LineWidth',2);
hold on;
plot3(X_Centroid_px,Y_Centroid_px,Z_Centroid_px,'MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','o');
view(-66,24);
rotate3d on;
hold on;