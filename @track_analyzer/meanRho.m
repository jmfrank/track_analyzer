%Calculates the mean number of cells over an experiment. Useful for
%short time lapses where cells aren't dividing a lot over the course of
%experiment. 
function [mean_rho_units,mean_rho] = meanRho(obj)
    


%Get exp info
exp_info = obj.exp_info;

%For each time point, get the total number of cell detections
for i = 1:length(obj.cells)
    %Number of cells
    n_cells = length( obj.cells{i} );

end

%Mean cells / time point
mean_rho = mean( n_cells );
%Mean normalized to area of image (cells / mm^2)
mean_rho_units = mean( n_cells ) ./ ( exp_info.pixel_size*exp_info.image_size(1)*exp_info.image_size(2))*10^6;

end




