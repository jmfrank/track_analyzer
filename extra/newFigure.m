
%Starts a new figure for nice plotting. 
function fig_h  = newFigure(i)

if(nargin < 1)
    fig_h = figure;
else
    fig_h = figure(i);
end

clf(fig_h);
set(fig_h,'color','w');
set(gca,'FontSize',20);
set(gca,'Linewidth',3);
box off;





end