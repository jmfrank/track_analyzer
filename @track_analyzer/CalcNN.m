
%Calculate the distance to the nearest neighbor. 
function [obj] = CalcNN(obj)
    


%Get exp info
exp_info = obj.exp_info;

%For each time point, estimate the local density of cells for each contour.
for i = 1:length(obj.cells)
    display(['Analyzing frame #',num2str(i)])
    %get x,y position of each cell
    dat = cell2mat(cellfun(@(x) mean(x), obj.cells{i},'UniformOutput',false)');
    D = pdist2(dat,dat);
    %Loop over cells in frame 'i'
    for j = 1:size(D,1)
        %Find the nearest neighbor 
        sel = D(j,:) > 0;
        nn(i,j)        = min( D(j,sel));
    end 
end

%Now go through the tracking data and add in the neighbor metric
for i = 1:length(obj.tracks)
    for j = 1:size(obj.tracks{i},1)
        t = obj.tracks{i}(j,1);
        c = obj.tracks{i}(j,2);
        obj.tracks{i}(j,7) = nn(t,c);
        
    end
end



end




