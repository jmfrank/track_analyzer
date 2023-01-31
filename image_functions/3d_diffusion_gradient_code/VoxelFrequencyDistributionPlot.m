%%  Written by Bhavna Rajasekaran, PhD
%%  Max Planck Institute for Physics of Complex Systems, Dresden, Germany
%%  Contact: bhavna@pks.mpg.de / bhavnarajasekaran@yahoo.com
%%  Reference: Object segmentation and ground truth in 3D embryonic imaging, Stem Plot Figure 2c,2d
%%  Refer to DS_GMM_Kmeans_3dnucleisegmentation.m

%%  Returns Frequency Distribution of Voxels plot for each x,y, z directions
%%  Requires 3D segmented image I with frequency list properties from DS_GMM_Kmeans_3dnucleisegmentation.m

function [Pxls_yrange, Pxls_xrange, Pxls_zrange, border]=VoxelFrequencyDistributionPlot(I,freqx,freqy,freqz)

border=2;
[nx,ny,nz]=size(I);
figure;

subplot(3,1,1); bar(freqx(2:end,1),freqx(2:end,2),0.5,'FaceColor',[0.419608 0.556863 0.137255],'EdgeColor','none'); %olivedrab
Pxls_yrange=max([1 (freqx(2,1)-border)]):min([(freqx(end,1)+border) ny]);
box off;
xlabel('X (pixels)','FontSize',12)
ylabel('Frequency','FontSize',12)
set(gca,'TickLength',[0 0],'XMinorTick','off','YMinorTick','off','YTick',[0 max(ylim)]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
hold on;

subplot(3,1,2); bar(freqy(2:end,1),freqy(2:end,2),0.5,'FaceColor',[0.576471 0.439216 0.858824],'EdgeColor','none'); %mediumpurple
Pxls_xrange=max([1 (freqy(2,1)-border)]):min([(freqy(end,1)+border) nx]);
box off;
xlabel('Y (pixels)','FontSize',12)
ylabel('Frequency','FontSize',12)
set(gca,'TickLength',[0 0],'XMinorTick','off','YMinorTick','off','YTick',[0 max(ylim)]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
hold on;

subplot(3,1,3);bar(freqz(2:end,1),freqz(2:end,2),0.5,'FaceColor',[1 0.54902 0],'EdgeColor','none');%orange
Pxls_zrange=max([1 (freqz(2,1)-border)]):min([(freqz(end,1)+border) nz]);
box off;
xlabel('Z (pixels)','FontSize',12)
ylabel('Frequency','FontSize',12)
set(gca,'TickLength',[0 0],'XMinorTick','off','YMinorTick','off','YTick',[0 max(ylim)]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
hold on;