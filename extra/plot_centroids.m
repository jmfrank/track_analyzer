%Plotting centroids. 

function plot_centroids(ids, ctrs)

hold on
for i = 1:length(ids)
    
    text(ctrs(ids(i),1),ctrs(ids(i),2),num2str(ids(i)),'color','g')

end
hold off


end