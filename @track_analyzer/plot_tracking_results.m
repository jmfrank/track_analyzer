%New plot tracking function for working with LSM images directly. 


function plot_tracking_results(obj,HA,step,params)



Z = obj.exp_info.z_planes;
T = obj.exp_info.t_frames;


%open the whole image
IMG = bfopen(obj.exp_info.img_file);

%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = bfGetReader(obj.exp_info.img_file);
series = 1;

%Get the image size of this series. 
size_x = reader.getSizeX;
size_y = reader.getSizeY;

%Before we plot, just do an average of the signal channel
mip = zeros(size_y,size_x,T);
for t = 1:T
    
    
 %Create empty image to fill up
    img = zeros(size_y,size_x,Z);
    
    %Get the bio-formats image index corresponding to this z-stack:
    for i = 1:Z
        this_plane = reader.getIndex(i-1,params.this_channel-1,t-1)+1;
        img(:,:,i) = IMG{series}{this_plane,1};
    end
    
    %Get the segment channel z-stack for this time:
    z_planes = 1:Z;
    
    %Max project
    mip(:,:,t) = max(img,[],3);
    
end

IDS = find(cellfun('length',obj.tracks) > params.min_length);

color_vec = linspecer(length(IDS));
%Now loop over time, but plot tracks with labels
figure(HA)

for t = 1:T
    
    
    imagesc(mip(:,:,t));
    colormap gray
    hold on
    
    %Need to figure out which tracks are present, then plot
    c=1;
    for j = 1:length(IDS)
        track = obj.tracks{IDS(j)};
        t_idx = find(track(:,1)==t);
        if(~isempty(t_idx))
            %centroids = track(1:t_idx,3:4);
            %C(i) = plot(centroids(:,1),centroids(:,2),'-*','color',color_vec(j,:),'Markersize',7,'linewidth',5);
            %label(C(i),num2str(id(j)),'color','w');
            %Find current cell mask
            cell_id = obj.tracks{IDS(j)}(t_idx,2);
            cell = obj.cells{ t }{cell_id};
            C(c) = plot(cell(:,1),cell(:,2),'color',color_vec(j,:),'linewidth',5);
            if(params.label_centroids)
                label(C(c), ['Track: ',num2str(IDS(j))]); %,' Centroid id: ',num2str(cell_id) ], 'color','r');
            end
            c=c+1;

        end
    %centroids=[];
    end
%     
%     for c = 1:length(cells)
%         
%         plot( cells{c}(:,1),cells{c}(:,2),'linewidth',2)
%     end
    title(['Frame: ',num2str(t)])
    hold off
    pause
    
    
    
end



end