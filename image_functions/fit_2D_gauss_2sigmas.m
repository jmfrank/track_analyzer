%Fit a 3D gaussian to 3D intensity data

%Input 3D image stack
%Centroid is [y,x] estimated center of blob
% y,x,z are the 


function fit = fit_2D_gauss_2sigmas( img, centroid, params)

%x is the fit to the 2D gaussian parameters are: 
%    x(1) = sigma-y
%    x(2) = sigma-x
%    x(3) = y-location
%    x(4) = x-location
%    x(5) = scaling factor
%    x(6) = background
%    a = y 
%    b = x 
%    c = z 

%Linearize stack
int_array = img(:);
%Get an array list for the coordinates
[I,J] = ind2sub(size(img),[1:length(int_array)]);
a = I';
b = J';

%Create 2D gaussian function to minimize:
fun = @(x) x(5).*exp( - (a - x(3)).^2./(2*x(1)^2) - (b - x(4)).^2./(2*x(2)^2)) + x(6) - int_array;

%Perform non-linear least squares:

    %Set up initial guess:
    x0(1) = params.xy_sigma_est;
    x0(2) = params.xy_sigma_est;
    x0(3) = centroid(2);
    x0(4) = centroid(1);
    x0(5) = params.int_scale;
    x0(6) = params.bg;
    
    %Bounds for optimizing
    lb = [0.2,0.2,1,1,10,10];
    ub = [2,2,size(img,1),size(img,2),5000,500];
    options = optimoptions('lsqnonlin','Display','off');

    [vals,resnorm] = lsqnonlin(fun,x0,lb,ub,options);
 
%Output structure
fit.sigma_y = vals(1);
fit.sigma_x = vals(2);
fit.y_pos   = vals(3);
fit.x_pos   = vals(4);
fit.int     = vals(5);
fit.bg      = vals(6);
fit.sigma_z = NaN;
fit.z_pos   = NaN;
fit.dist_2_periphery = NaN;
fit.cell_id = NaN;
fit.resnorm = resnorm;
    %We also want the total intensity above background from the raw data
    sel = int_array > fit.bg;
    fit.raw_sum = sum( int_array(sel) );
    

end




