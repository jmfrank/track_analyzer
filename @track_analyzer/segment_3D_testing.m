% 3D Segmentation of cell nuclei. Based on: 'Object segmentation and ground
% truth in 3D embryonic imaging'. 

% Specify segmentation parameters in params, and function options in step.
% Written for single time-point images. <<< will write for looping over
% time later...

function segment_3D_testing(img, params, step)

%Get debug
debug=step.debug;

if(debug)
    params
end

%% Pre-processing using bio-formats. 

%Get the image size of this series.    
T = 1;
size_x = size(img,2);
size_y = size(img,1);


%% Generate im_info structure
    ZSlicesinStack = size(img,3);
    xscale         = 1;
    yscale         = 1;
    zscale         = 1;
    image_bits     = 16;
    
    
%% Segmentation parameters
    AbsMinVol = params.AbsMinVol;
    AbsMaxVol = params.AbsMaxVol; %Maximum voxel per single nucleus

    
    
%% Segmentation steps
    display_plots = step.display_plots;
    
    
%% Set Segmentation Algorithm Parameter values
    %%Derivatives Sum (DS) Algorithm
    alpha = params.alpha;
    beta = params.beta;
    gamma = params.gamma;
    epsilon =params.epsilon;
    delta = params.delta;
    sigmagradient=params.sigmagradient;

    %%Noise filter parameters
    h = params.noise_filter; 
    w1=params.median_filter(1);w2=params.median_filter(2); %%Median filter

    %%Perona-Malik based non-linear isotropic diffusion filter
    diffuse_iterations=params.diffuse_iterations;
    kappa1=params.kappa1;
    kappa2=params.kappa2;
    option=params.option;
    
    %Thresholding
    thrshlevel = params.thrshlevel;
    % nconn value for neighborhood calculations
    nconn_BW = params.nconn_BW;
    nconn_BW2= params.nconn_BW2;
    
    boundary_param = params.boundary_param;

%% Loop over images
for t=1:T
    
    %Get image planes for this time point. 
    I = zeros(size_y,size_x,ZSlicesinStack);
    %Empty image matrices.
    G = I; Diff_im = I; Fim = I; Mag_fim = I; Mag_Lim = I; Det_hessian = I;
    Sum = I; Multi = I; J = I;
     
    tic
    %Get the bio-formats image index corresponding to this z-stack:
    for i = 1:ZSlicesinStack
        this_plane     = reader.getIndex(i-1,params.seg_channel-1,t-1)+1;
        this_plane_img = bfGetPlane(reader,this_plane);
        I(:,:,i)     = this_plane_img;
    end
    
    toc
    
    
 
    for zslice= 1:ZSlicesinStack %% iterate through number of slices
        %
        %%use imshow to plot any of the Steps from 1-5. eg.figure; imshow((I(:,:,zslice))
        
        %%Step1: Denoising Filters
        %G(:,:,zslice) = imfilter(I(:,:,zslice), h,'replicate');
        %%Optional filters median or deconvlucy
        G(:,:,zslice) = medfilt2(I(:,:,zslice), [w1 w2]);
        %G(:,:,zslice) = deconvlucy(G(:,:,zslice),h); %%use the same 'h' as for Gaussian or define new variable
        
        %%Step2: Perona & Malik non-linear isotropic diffusion filter, refer to
        %%diffusioncode.m for details, here alter only diffuse_iterations & kappa1       
        Diff_im(:,:,zslice) = diffusioncode(G(:,:,zslice), diffuse_iterations, 0.1429, kappa1, kappa2, option);
        
        %%Step3a: Gauss gradient, 1st derivative of the image, refer to gaussgradient.m for details
        Fim(:,:,zslice)=mat2gray(Diff_im(:,:,zslice));
        [imx,imy]=gaussgradient(Fim(:,:,zslice),sigmagradient);
        Mag_fim(:,:,zslice) = hypot(imx,imy);
        
        %%Step3b: Laplacian
        [L_imxx, L_imxy] = gaussgradient(imx,sigmagradient);
        [L_imyx, L_imyy] = gaussgradient(imy,sigmagradient);
        Mag_Lim(:,:,zslice) = L_imxx+L_imyy;  %Laplacian (trace of the 2D matrix)
        Mag_Lim(:,:,zslice)=Mag_Lim(:,:,zslice).*(Mag_Lim(:,:,zslice)>0);
        
        %%Step3c: Hessian Determniant
        Det_hessian(:,:,zslice)=L_imxx.*L_imyy-L_imxy.*L_imyx;
        Det_hessian(:,:,zslice)=Det_hessian(:,:,zslice).*(Det_hessian(:,:,zslice)<0);
        %%Step4: Masking function using tanh on the Summed Derivatives from Step3
        Sum(:,:,zslice)=1-(alpha*Mag_fim(:,:,zslice)+beta*Mag_Lim(:,:,zslice)+epsilon*(abs(Det_hessian(:,:,zslice))));
        X=gamma-((alpha*Mag_fim(:,:,zslice) + beta*Mag_Lim(:,:,zslice) + epsilon*abs(Det_hessian(:,:,zslice))) /delta);
        Multi(:,:,zslice)=0.5*(tanh(X)+1);% masking function
        J(:,:,zslice)=(double(G(:,:,zslice)).*Multi(:,:,zslice)); % masked image applied on smoothened image
        

    end
    
    if(debug)
        imshow3D(J)
        pause
    end
    
    %% Global Otsu's method. Image Thresholding using Otsu's method, 8-bit or 16-bit
      
    BW = zeros(size(I));
    if (image_bits==8)
        %thrshlevel = otsuthresh( uint8(J(:)));
        BW(:,:,zslice)=im2bw(uint8(J(:,:,zslice)),thrshlevel);
    elseif(image_bits==16)
        %thrshlevel = otsuthresh( uint16(J(:)));
        for i = 1:1:ZSlicesinStack
            %Threshold 
            BW(:,:,i) = im2bw(uint16(J(:,:,i)),thrshlevel);
            %Also fill holes laterally. Vertical holes can be
            %problematic if we just to imfill with 3D stack
            BW(:,:,i) = imfill(BW(:,:,i),'holes');
        end

    end
    

    %% Filling image. Let's try bounary.     
    CC = bwconncomp(logical(BW), nconn_BW);
    stats = regionprops(CC,'Area','PixelList','PixelIdxList');
    
    if(debug)
        imshow3D(BW)
        pause
    end
    
    %% Fill in gaps using 'boundary' and 'voxelization' using sub-function vox_this_frame. 
    
    BW2 = vox_this_frame(BW, stats, params);
    disp(['Finished voxing.'])
    %% Plotting result for debugging
    if(debug)
        figure(1); clf
        imshow3D(BW2)
%         CC = bwconncomp(logical(BW2), 6);
%         stats = regionprops(CC,'Area','PixelList');
%         for i = 1:length(stats)
%             ps = stats(i).PixelList;
%             ps = ps(:,[2,1,3]);
%             shp = alphaShape(ps);
%             plot(shp)
%             pause
%         end
    end
        
    
    %% Step6: Connect similar pixels of zslices to get a 3D binary stack
    
    %%Step7: Computes three properties for each 3D segmented object
    Con=bwconncomp(logical(BW2),nconn_BW2);
    stats = regionprops(Con,'Centroid','Area','PixelList','PixelIdxList');
    numobj=numel(stats); %total number of segmented objects
    
    %%Step8: Image stack dimensions & image scaling in mirons
    [nx,ny,nz]=size(I);
    x=(1:nx)*xscale;
    y=(1:ny)*yscale;
    z=(1:nz)*zscale;
    %%Basic DS Algorithm Segmentation acheived
    
%     %% Objects Volume analysis to identify outliers (potentially fused objects)    
%     volume_vector=0;
%     counter=1;
%     for s=1:length(stats)
%         stats_volume_vector(s)=stats(s).Area;
%         if (stats(s).Area>AbsMinVol && stats(s).Area<AbsMaxVol)
%             volume_vector(counter)=stats(s).Area; % an array of volumes of objects
%             counter=counter+1;
%         end
%     end
%     
%     
%     %%Determine fused object volume threshold & seggregate outliers
%     outliercut=0.1;
%     outliers=volume_vector-mean(volume_vector)>(outliercut*std(volume_vector)); % determine fused objects
%     outlier_index=find(outliers);
%     outlier_volume_vector=volume_vector(outlier_index);
%     more_outliers=find(stats_volume_vector>=AbsMaxVol); %more fused nuclei
%     final_outlier_volume=horzcat(outlier_volume_vector,stats_volume_vector(more_outliers)); % outliers for post-processing
%     
%     liar_list=find(ismember(stats_volume_vector,final_outlier_volume));%find the index of outliers in the original stats based on volume criteria
%     mean_volume=mean(volume_vector);
%     std_volume=std(volume_vector);
%     %%Emperically determined oultier volume cut -off using mean_volume & standard deviation_volume
%     outliers_volume_cutoff=mean_volume+(outliercut*std_volume);
%     
%     %%About noisy voxels below 10 voxels that are discarded
%     noisy_outliers=find(stats_volume_vector<=AbsMinVol);
%     if( ~isempty(noisy_outliers) )
%         noisy_list=find(ismember(stats_volume_vector,stats_volume_vector(noisy_outliers)))
% 
%         for l=1:length(noisy_list)
%             noisy_volume(l)=stats(noisy_list(l)).Area;
%             noisy_volume_microns_cube(l)=noisy_volume(l)*xscale*yscale*zscale;
%         end
%         maximum_noisy_volume_microns_cube=max(noisy_volume_microns_cube);
%     end
%     
% 
%     %% Seggregate correctly segmented objects by the DS ALgorithm (voxels>noise(10) & voxels<Maximum_allowed_volume_size)
%     Maximum_allowed_volume_size = min(final_outlier_volume);
%     Minimum_allowed_volume_size = AbsMinVol;
%     cn=1;
%     for correctvol=1:length(stats)
%         if (stats(correctvol).Area < Maximum_allowed_volume_size && (stats(correctvol).Area > Minimum_allowed_volume_size))
%             Centroid_px(cn,1)=stats(correctvol).Centroid(1)*xscale;
%             Centroid_px(cn,2)=stats(correctvol).Centroid(2)*yscale;
%             Centroid_px(cn,3)=stats(correctvol).Centroid(3)*zscale;
%             Pix_lis{cn}(:,1)=stats(correctvol).PixelList(:,1)*xscale;
%             Pix_lis{cn}(:,2)=stats(correctvol).PixelList(:,2)*yscale;
%             Pix_lis{cn}(:,3)=stats(correctvol).PixelList(:,3)*zscale;
%             cn=cn+1;
%         end
%     end
%     
    %% Allocate Correctly Segmented Objects by the DS algorithm
%     for r=1:size(Centroid_px,1)
%         %%Correctly segmented DS Algorithm
%         Final_cen_px(r)=struct('Centroid_xyz',Centroid_px(r,1:3),'Voxellistvalues',Pix_lis{r});
%         %%Initialize Final_cen_px_kmeans(DS+Kmeans), assign DS results
%         Final_cen_px_kmeans(r)=struct('Centroid_xyz',Centroid_px(r,1:3),'Voxellistvalues',Pix_lis{r});
%         %%Initialize Final_cen_px_gmm (DS+GMM), assign DS results
%         Final_cen_px_gmm(r)=struct('Centroid_xyz',Centroid_px(r,1:3),'Voxellistvalues',Pix_lis{r});
%     end
%    disp('The number of objects correctly segmented without post processing by DS algorithm: ');
%    disp(length(Centroid_px));
    
    %% Post-Processing of Fused Objects with GMM & Kmeans
     %%By setting to zero, you can switch of post-processing of fused nuclei
    if step.fused_nuclei
        
        %%setting scales to find next peak for fused nuclei based on
        %%image resolution in x,y,z (assessed from micrsocope resolution)
        %%used for voxel frequency distribution
        scale_peak_distx=round(2.768/xscale);
        scale_peak_disty=round(2.768/yscale);
        scale_peak_distz=round(2/zscale);
        
        noteindex_GMM =[];%%initialize a variable to track indices of objects that remain fused after GMM,
        %%in case many are fused after GMM, then change parametrs of noise
        %%filters & DS algorithm to improve segmentation accuarcy
        
        %%post-process each fused object from liar_list variable
        for l=1:length(liar_list)
            
            liar_Area=stats(liar_list(l)).Area;
            Pxlist=stats(liar_list(l)).PixelList;
            Px=Pxlist;
            
            %%micron scaling of pixels
            Px(:,1)=Px(:,1)*xscale;
            Px(:,2)=Px(:,2)*yscale;
            Px(:,3)=Px(:,3)*zscale;
            
            %%micron scaling of Centroids
            X_Centroid_px=(stats(liar_list(l)).Centroid(1))*xscale;
            Y_Centroid_px=(stats(liar_list(l)).Centroid(2))*yscale;
            Z_Centroid_px=(stats(liar_list(l)).Centroid(3))*zscale;
            
            Pxlist(end+1,1:3)=[-1 -1 -1]; % to avoid index out of bound error
            
            %%Determine frequency distribution of voxels in x, y, z
            %%directions.Use inter-peak distances within voxel frequency
            %%distribution to determine number of fused nuclei
            
            %%X-Direction frequency distribution
            freqx=tabulate(Pxlist(:,1));
            if (length(freqx(:,1))>10)
                if (length(freqx(2:end,2))>5), peak_distx=scale_peak_distx;else peak_distx=1;end
                if peak_distx==1
                    x_peaks=0;
                else
                    x_peaks=numel(findpeaks(freqx(2:end,2),'minpeakdistance',peak_distx));
                end
            else
                x_peaks=0;
            end
            
            %%Y-Direction frequency distribution
            freqy=tabulate(Pxlist(:,2));
            if (length(freqy(:,1))>10)
                if (length(freqy(2:end,2))>5), peak_disty=scale_peak_disty;else peak_disty=1;end
                if peak_disty==1
                    y_peaks=0;
                else
                    y_peaks=numel(findpeaks(freqy(2:end,2),'minpeakdistance',peak_disty));
                end
            else
                y_peaks=0;
            end
            
            %%Z-Direction frequency distribution
            freqz=tabulate(Pxlist(:,3));
            if (length(freqz(:,1))>5)
                peak_distz=scale_peak_distz;
                if (length(freqz(2:end,2))>5),z_peaks=numel(findpeaks(freqz(2:end,2),'minpeakdistance',peak_distz)); else z_peaks=0;end
            else
                z_peaks=0;
            end
            
            %% Plot VoxelFrequencyDistributionPlot
            if display_plots
                VoxelFrequencyDistributionPlot(I,freqx,freqy,freqz); %% Plot of frequency distribution of voxels here
                subplot(3,1,1);title('Voxel frequency distribution for finding peaks for GMM & Kmeans','FontSize',10);
            end
            %%
            
            %% Render nuclei in 3D for the respective nuclei under post-processing
            if display_plots
                RenderNucleiDS(G,x,y,z,freqx,freqy,freqz,X_Centroid_px,Y_Centroid_px,Z_Centroid_px);
            end
            
            %%determine possible number of fused nuclei (maximum of peaks obtained from x,y,z directions)
            nuclei_cluster=max([x_peaks y_peaks z_peaks]);
            
            %% GMM Post-Processing/Local Clustering of fused Candidate Nuclei with GMM
            %By setting to zero, you can switch off GMM
            if step.Active_GMM
                
                if (length(freqx) > 5 && length(freqy) > 5 && length(freqz) > 5)
                    nuclei_cluster_add=nuclei_cluster+4; %% start with more clusters than findpeaks gives to intialze gmm
                    AIC = zeros(1,nuclei_cluster_add);
                    options = statset('MaxIter',800, 'TolFun',1e-8,'Display','Off');
                    
                    for k = 1:nuclei_cluster_add
                        obj{k} = gmdistribution.fit(Px,k,'Start','randSample','Replicates',5,'Regularize',0.2,'SharedCov',true,'CovType','diagonal','Options',options);
                        AIC(k)= obj{k}.AIC;
                    end
                    [minAIC,numComponents] = min(AIC); %%AIC dtermines best fitting model to the voxel data
                    model = obj{numComponents}; %% Choose the GMM model determined by AIC
                    
                    %%Check if AIC determined GMM has not over-split nuclei
                    if (numel(model.mu(:,1))>1)
                        D = pdist(model.mu);
                        distancecut=3.5; % use a cut-off to check from over-splitting nuclei whose centroids are <3.5 micrometer distance
                        if (length(find(D<=distancecut)))>=1
                            outlier_centroids=numel(model.mu(:,1))-length(find(D<=distancecut)); %% find no. of centroids that are less than 3.5 microns distance
                            if (outlier_centroids <= numComponents && outlier_centroids>0),numComponents=outlier_centroids;end
                            if outlier_centroids==0;numComponents=1;end
                        end
                    end
                    
                    model = obj{numComponents};
                    for comp=1:numComponents
                        [idxgm,nlogl,Prob,logpdf,P] = cluster(obj{comp},Px); %%find index of which pixels belong to which cluster of nuclei
                    end
                    
                    %%Assign New Centroids & Pixels and concatenate results
                    new_cluster_cat=[];
                    for clusterform=1:max(idxgm)
                        new_cluster = Px(idxgm == clusterform,:);
                        new_cluster_cat=cat(1,new_cluster_cat, new_cluster);
                        new_centroids=model.mu(clusterform,:);
                        gmm_cluster(clusterform)=struct('Centroid_xyz',new_centroids(:,1:3),'Voxellistvalues',new_cluster(:,1:3));
                        if (length(new_cluster)>400) %check for voxels>400 after GMM, or your choice that are surely fused nuclei, if noteindex_GMM=1, good result
                            noteindex_GMM=cat(2,noteindex_GMM,l);
                        end
                    end
                    Final_cen_px_gmm=cat(2,Final_cen_px_gmm,gmm_cluster);
                    
                    if display_plots
                        GMMCentroids(model);
                    end
                    clear gmm_cluster new_cluster new_centroids model numComponents obj AIC outlier_centroids idxgm;
                else
                    %%Not process through GMM if the nucleus frequency list<5 in length
                    C(:,1)= X_Centroid_px;
                    C(:,2)= Y_Centroid_px;
                    C(:,3)= Z_Centroid_px;
                    clustered_once(1)=struct('Centroid_xyz',C(1,1:3),'Voxellistvalues',Px(:,1:3));
                    Final_cen_px_gmm=cat(2,Final_cen_px_gmm,clustered_once);
                    clear clustered_once;
                end
            end % end GMM clsuering
            
            %% Kmeans Post-Processing/Local Clustering of fused Candidate Nuclei using Kmeans
            %By setting to zero, you can switch off Kmeans
            IDX=0;% for Kmeans to check number of clusters
            if step.Active_Kmeans
                
                if (nuclei_cluster>=1)
                    [IDX,C] = kmeans(Px,nuclei_cluster,'distance','sqEuclidean','onlinephase','on','emptyaction','drop','replicates',5);
                    
                    if  (max(IDX)==1) %%
                        clustered_once(1)=struct('Centroid_xyz',C(1,1:3),'Voxellistvalues',Px(:,1:3));
                        Final_cen_px_kmeans=cat(2,Final_cen_px_kmeans,clustered_once);
                        clear clustered_once;
                    end
                else
                    C(:,1)= X_Centroid_px;
                    C(:,2)= Y_Centroid_px;
                    C(:,3)= Z_Centroid_px;
                    clustered_once(1)=struct('Centroid_xyz',C(1,1:3),'Voxellistvalues',Px(:,1:3));
                    Final_cen_px_kmeans=cat(2,Final_cen_px_kmeans,clustered_once);
                    clear clustered_once;
                end
                
                if display_plots
                    KmeansCentroids(C);
                    hold off;
                end
                
                if IDX
                    %%Prepare Kmeans for second iteration
                    if (max(IDX)>1)
                        for ind=1:max(IDX)
                            clear row_num NewPxlist New_Px;
                            row_num=find(IDX==ind);
                            
                            NewPxlist=Pxlist(row_num,1:3); %% in Pixel units
                            New_Px=Px(row_num,1:3);  %% scaled in microns
                            NewPxlist(end+1,1:3)=[-1 -1 -1];
                            nborder=3;
                            
                            %%Calculate Frequency distribtions of voxels in each
                            %%direction for the new voxel list obtained from Kmeans
                            %%previous iteration
                            
                            %X-direction voxel frequency distribution
                            nfreqx=tabulate(NewPxlist(:,1));
                            if (length(nfreqx(:,1))>10)
                                if (length(nfreqx(2:end,2))>5), npeak_distx=scale_peak_distx;else npeak_distx=1;end
                                if npeak_distx==1
                                    nx_peaks=0;
                                else
                                    nx_peaks=numel(findpeaks(nfreqx(2:end,2),'minpeakdistance',npeak_distx));
                                end
                            else
                                nx_peaks=0;
                            end
                            
                            %Y-direction voxel frequency distribution
                            nfreqy=tabulate(NewPxlist(:,2));
                            if (length(nfreqy(:,1))>10)
                                if (length(nfreqy(2:end,2))>5), npeak_disty=scale_peak_disty;else npeak_disty=1;end
                                if npeak_disty==1
                                    ny_peaks=0;
                                else
                                    ny_peaks=numel(findpeaks(nfreqy(2:end,2),'minpeakdistance',npeak_disty));
                                end
                            else
                                ny_peaks=0;
                            end
                            
                            %Z-direction voxel frequency distribution
                            nfreqz=tabulate(NewPxlist(:,3));
                            if (length(nfreqz(:,1))>5)
                                peak_distz=scale_peak_distz;
                                if (length(nfreqz(2:end,2))>5), nz_peaks=numel(findpeaks(nfreqz(2:end,2),'minpeakdistance',peak_distz)); else nz_peaks=0;end
                                
                            else
                                nz_peaks=0;
                            end
                            %% Plot VoxelFrequencyDistributionPlot for K-means 2nd iteration
                            %if display_plots
                            %   VoxelFrequencyDistributionPlot(I,nfreqx,nfreqy,nfreqz) %% Plot of frequency distribution of voxels here
                            %  subplot(3,1,1);title('Voxel frequency distribution for Kmeans 2nd iteration','FontSize',10);
                            %end

                            new_cluster=max([nx_peaks ny_peaks nz_peaks]);%%Find maximum of peaks
                            
                            if (new_cluster>1)
                                [new_ID,More_C] = kmeans(New_Px,new_cluster,'distance','sqEuclidean','onlinephase','on','emptyaction','drop','replicates',5);
                                
                                if More_C, C(ind,1:3)=[0 0 0]; end
                                
                                if (max(new_ID)==1)
                                    clustered_one(1)=struct('Centroid_xyz',More_C(1,1:3),'Voxellistvalues',New_Px(:,1:3));
                                    Final_cen_px_kmeans=cat(2,Final_cen_px_kmeans,clustered_one);
                                    clear clustered_one;
                                elseif (max(new_ID)>1)
                                    for clust_id2=1:max(new_ID)
                                        cluster_index2=find(new_ID==clust_id2);
                                        clear Reassign_px2;
                                        Reassign_px2{clust_id2}=New_Px(cluster_index2,1:3);
                                        clustered_onceagain(clust_id2)=struct('Centroid_xyz',More_C(clust_id2,1:3),'Voxellistvalues',Reassign_px2{clust_id2});
                                    end
                                    Final_cen_px_kmeans=cat(2,Final_cen_px_kmeans,clustered_onceagain);
                                    clear clustered_onceagain;
                                end
                            else
                                clear More_C;
                                More_C(:,1)= C(ind,1);
                                More_C(:,2)= C(ind,2);
                                More_C(:,3)= C(ind,3);
                                clear not_clustered;
                                not_clustered(1)=struct('Centroid_xyz',C(ind,1:3),'Voxellistvalues',New_Px(:,1:3));
                                Final_cen_px_kmeans=cat(2,Final_cen_px_kmeans,not_clustered);
                            end
                        end
                    end
                else
                end
            end % end Kmeans clsuering
            
            if display_plots dump=input(['y']);close all; end
            
        end
    end % end of fused nuclei processing

%     %% Segmentation Result from DS Algorithm
%     nuclei_centroid_voxels_microns{t}=Final_cen_px;
%     nuclei_centroid_microns_px_voxels{t}=nuclei_centroid_voxels_microns{t};
%     for tot_n=1:length(nuclei_centroid_voxels_microns{t})   %% for all nuclei in a stack, centroids & voxels both in microns
%         for px_len=1:length(nuclei_centroid_voxels_microns{t}(tot_n).Voxellistvalues(:,1))   %% for each object in that stack, centroid in microns, voxels in voxel units
%             nuclei_centroid_microns_px_voxels{t}(tot_n).Voxellistvalues(:,1)=round(nuclei_centroid_voxels_microns{t}(tot_n).Voxellistvalues(:,1)/xscale);
%             nuclei_centroid_microns_px_voxels{t}(tot_n).Voxellistvalues(:,2)=round(nuclei_centroid_voxels_microns{t}(tot_n).Voxellistvalues(:,2)/yscale);
%             nuclei_centroid_microns_px_voxels{t}(tot_n).Voxellistvalues(:,3)=round(nuclei_centroid_voxels_microns{t}(tot_n).Voxellistvalues(:,3)/zscale);
%         end
%     end
%     
%     %% Segmentation Result from DS + GMM Algorithm
%     nuclei_centroid_voxels_microns_gmm{t}=Final_cen_px_gmm;
%     nuclei_centroid_microns_px_voxels_gmm{t}=nuclei_centroid_voxels_microns_gmm{t};
%     for tot_n=1:length(nuclei_centroid_voxels_microns_gmm{t})   %% for all nuclei in a stack, centroids & voxels both in microns
%         for px_len=1:length(nuclei_centroid_voxels_microns_gmm{t}(tot_n).Voxellistvalues(:,1))   %% for each object in that stack, centroid in microns, voxels in voxel units
%             nuclei_centroid_microns_px_voxels_gmm{t}(tot_n).Voxellistvalues(:,1)=round(nuclei_centroid_voxels_microns_gmm{t}(tot_n).Voxellistvalues(:,1)/xscale);
%             nuclei_centroid_microns_px_voxels_gmm{t}(tot_n).Voxellistvalues(:,2)=round(nuclei_centroid_voxels_microns_gmm{t}(tot_n).Voxellistvalues(:,2)/yscale);
%             nuclei_centroid_microns_px_voxels_gmm{t}(tot_n).Voxellistvalues(:,3)=round(nuclei_centroid_voxels_microns_gmm{t}(tot_n).Voxellistvalues(:,3)/zscale);
%         end
%     end
%     if (~isempty(nuclei_centroid_microns_px_voxels_gmm{1}))
%         for m=1:length(nuclei_centroid_microns_px_voxels_gmm{1})
%             segmentedvolumes_gmm(m)=length(nuclei_centroid_microns_px_voxels_gmm{1}(m).Voxellistvalues);
%         end
%     else
%         for m=1:length(nuclei_centroid_microns_px_voxels_gmm{2})
%             segmentedvolumes_gmm(m)=length(nuclei_centroid_microns_px_voxels_gmm{2}(m).Voxellistvalues);
%         end
%     end
%     
%     %% Segmentation Result from DS + Kmeans Algorithm
%     nuclei_centroid_voxels_microns_kmeans{t}=Final_cen_px_kmeans;
%     nuclei_centroid_microns_px_voxels_kmeans{t}=nuclei_centroid_voxels_microns_kmeans{t};
%     for tot_n=1:length(nuclei_centroid_voxels_microns_kmeans{t})   %% for all nuclei in a stack, centroids & voxels both in microns
%         for px_len=1:length(nuclei_centroid_voxels_microns_kmeans{t}(tot_n).Voxellistvalues(:,1))   %% for each object in that stack, centroid in microns, voxels in voxel units
%             nuclei_centroid_microns_px_voxels_kmeans{t}(tot_n).Voxellistvalues(:,1)=round(nuclei_centroid_voxels_microns_kmeans{t}(tot_n).Voxellistvalues(:,1)/xscale);
%             nuclei_centroid_microns_px_voxels_kmeans{t}(tot_n).Voxellistvalues(:,2)=round(nuclei_centroid_voxels_microns_kmeans{t}(tot_n).Voxellistvalues(:,2)/yscale);
%             nuclei_centroid_microns_px_voxels_kmeans{t}(tot_n).Voxellistvalues(:,3)=round(nuclei_centroid_voxels_microns_kmeans{t}(tot_n).Voxellistvalues(:,3)/zscale);
%         end
%     end
%     if (~isempty(nuclei_centroid_microns_px_voxels_kmeans{1}))
%         for m=1:length(nuclei_centroid_microns_px_voxels_kmeans{1})
%             segmentedvolumes_kmeans(m)=length(nuclei_centroid_microns_px_voxels_kmeans{1}(m).Voxellistvalues);
%         end
%     else
%         for m=1:length(nuclei_centroid_microns_px_voxels_kmeans{2})
%             segmentedvolumes_kmeans(m)=length(nuclei_centroid_microns_px_voxels_kmeans{2}(m).Voxellistvalues);
%         end
%     end
%     
%     %% Plot SegmentedVolumesStemPlot
%     if display_plots
%         SegmentedVolumesStemPlot(mean_volume,std_volume,outliercut,volume_vector,segmentedvolumes_kmeans,segmentedvolumes_gmm)%%Stem Plot of Segmented Volumes
%     end
    
    
    
    %% Create frame_obj and save. 
    clear frame_obj
    %Add pixel list/centroid for each region to frame_obj
    counter = 1;
    for i = 1:length(stats)
        if (length(stats(i).PixelIdxList) >= AbsMinVol)
            frame_obj.PixelIdxList{counter} = stats(i).PixelIdxList;
            frame_obj.centroids{counter}    = stats(i).Centroid;
            counter = counter + 1;
        end
    end
    
    %Add final binarized image to frame_obj for save keeping
    frame_obj.BW = BW2;
    
%     %Lastly, keep the original output variables. Might need them later.
%         %Structured output from 'Object segmentation and ground. 
%         frame_obj.nuclei_centroid_voxels_microns = nuclei_centroid_voxels_microns;
%         frame_obj.nuclei_centroid_voxels_microns_gmm = nuclei_centroid_voxels_microns_gmm;
%         frame_obj.nuclei_centroid_voxels_microns_kmeans = nuclei_centroid_voxels_microns_kmeans;
%         
    %Save frame_obj as done before. 
    fname = ['frame_',sprintf('%04d',t),'.mat'];
    if(~exist([obj.exp_info.nuc_seg_dir],'dir'))
        mkdir(obj.exp_info.nuc_seg_dir)
    end
    save([obj.exp_info.nuc_seg_dir,fname],'frame_obj')
    disp(['Writing: ',fname] )

    
        
    
end

    

end




