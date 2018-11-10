% Calculate mean yap signal (could be arbitrary signal based on image).
% Uses the roi2poly tool. Maybe faster?
%Now computes cyto to nuclear ratio by taking a mean of a thin ring outside
%of the contour

function obj = calc_nuc_area(obj)


length(obj.tracks);
length(obj.img_files);


%Loop over tracks
for i = 1:length(obj.tracks)
    obj.nuc_area_trace{i} = [];
    
    %Loop over time points
    for j = 1:size(obj.tracks{i},1)
        
        %Time
        t = obj.tracks{i}(j,1);
        %cell ID
        id = obj.tracks{i}(j,2);
        
        %Find the cell
        x = obj.cells{t}{id}(:,1);
        y = obj.cells{t}{id}(:,2);
        
        %Calculate area
        A = polyarea(x,y);

        %Add to object
        obj.nuc_area_trace{i}(t) = A;        
        
    end
    
end


obj.nuc_area_calc = true;

end


% 
%         
% 
% for i = [1:length(obj.img_files)]
%     display(['Analyzing mean signal of frame: ',num2str(i)])
%     i
%     
%     for j = [1:length(obj.tracks)]
%         j
%         %if track exists at frame i
%         times = obj.tracks{j}(:,1)
%         
%         if( sum( times == i) == 1 ) 
%             time_idx = find(times == i);
%             cell_id = obj.tracks{j}(time_idx,2);
%             
%             %Contour
%             x     = obj.cells{i}{cell_id}(:,1);
%             y     = obj.cells{i}{cell_id}(:,2);
% 
%             %Calculate area
%             A = polyarea(x,y);
% 
%             %Add to object
%             A
%             obj.nuc_area_trace{j}(time_idx) = A;
%         end
%                 pause
% 
%     end
%     
% end
