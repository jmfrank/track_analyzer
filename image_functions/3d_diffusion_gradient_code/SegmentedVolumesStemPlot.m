%%  Written by Bhavna Rajasekaran, PhD
%%  Max Planck Institute for Physics of Complex Systems, Dresden, Germany
%%  Contact: bhavna@pks.mpg.de / bhavnarajasekaran@yahoo.com
%%  Reference: Object segmentation and ground truth in 3D embryonic imaging, Stem Plot Figure 2a
%%  Refer to DS_GMM_Kmeans_3dnucleisegmentation.m

function SegmentedVolumesStemPlot(mean_volume,std_volume,outliercut,volume_vector,segmentedvolumes_kmeans,segmentedvolumes_gmm)

volmax=1000;
nbins=40;
binsize=volmax/nbins;
bingrid=linspace(binsize/2,(volmax-(binsize/2)),nbins);
offset=binsize/6;

[nsegemented,xoutsegemented] = hist(volume_vector,bingrid); %%DS algorithm
[nkmeans,xoutkmeans] = hist(segmentedvolumes_kmeans,bingrid); %%DS+kmeans
[ngmm,xoutgmm] = hist(segmentedvolumes_gmm,bingrid); %%DS+GMM

figure;
h1=stem(xoutsegemented,nsegemented,'fill','Color',[0.651  0.361  0.643],'LineWidth',1,'MarkerEdgeColor',[0.651  0.361  0.643],'MarkerFaceColor',[0.651  0.361  0.643],'MarkerSize',2,'Marker','*');
hold on;
h3=stem(xoutkmeans-offset,nkmeans,'Color',[0.161  0.337  0.655],'LineWidth',1,'MarkerEdgeColor',[0.161  0.337  0.655],'MarkerFaceColor',[0.161  0.337  0.655],'MarkerSize',2,'Marker','o');
hold on;
h4=stem(xoutgmm+offset,ngmm,'Color',[0.435  0.745  0.267],'LineWidth',1,'MarkerEdgeColor',[0.435  0.745  0.267],'MarkerFaceColor',[0.435  0.745  0.267],'MarkerSize',2,'Marker','^');
hold on;
ylim('auto')
set(gca,'yscale','log');
hold on;
yaxis=0.1:1:max(ylim);
hold on;
x3=(mean_volume+outliercut*std_volume)*ones(1,max(ylim));
hold on;
plot(x3,yaxis,'--k','LineWidth',1);
hold on;
xlabel('Volumes of Objects (Voxels) ','FontSize',14)
ylabel('Number of Objects','FontSize',14)
set(gca,'FontSize',14);
set(gca,'TickLength',[0 0]);
set(gca,'LineWidth',1);
box off;
%title(['Segmented Volumes', ', Timepoint = ', num2str(time)],'FontSize',20);
leghdl=legend('DS algorithm','DS+Kmeans','DS+GMM', 'Fused object volume threshold');
legend('boxoff')
set(leghdl,'FontSize',14);
hold off;