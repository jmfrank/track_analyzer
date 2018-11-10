%Fit a 3D gaussian to 3D intensity data, but fix the xy_sigma and the
%z_sigma to the initial guess values [params.xy_sigma_est,
%params.z_sigma_est]. 

%Input 3D image stack
%Centroid is [y,x] estimated center of blob
% y,x,z are the 

function fit = fit_3D_gauss_fixed_sigma( stack, centroid, params)

%x is the fit to the 3D gaussian parameters are: 
%    x(1) = y-location
%    x(2) = x-location
%    x(3) = z-location
%    x(4) = scaling factor
%    x(5) = background
%    a = y 
%    b = x 
%    c = z 

sigma_xy = params.xy_sigma_est;
sigma_z  = params.z_sigma_est;

%Linearize stack
int_array = stack(:);
%Get an array list for the coordinates
[I,J,K] = ind2sub(size(stack),[1:length(int_array)]);
a = I';
b = J';
c = K';

%Create 3D gaussian function to minimize:
fun = @(x) x(4).*exp( - (a - x(1)).^2/(2*sigma_xy^2) - (b - x(2)).^2/(2*sigma_xy^2) - (c - x(3)).^2/(2*sigma_z^2)) + x(5) - int_array;

%Perform non-linear least squares:

    %Set up initial guess:
    x0(1) = centroid(2);
    x0(2) = centroid(1);
    x0(3) = centroid(3); 
    x0(4) = params.int_scale;
    x0(5) = params.bg;
        
    %Lower bound
    lb = [1,1,1,1,0];
    %Upper bound
    ub = [size(stack,1),size(stack,2),size(stack,3),65000,20000];
    options = optimoptions('lsqnonlin','Display','off');

    [vals, resnorm] = lsqnonlin(fun,x0,lb,ub,options);

%Correlation-coefficient
real_fun = @(x) x(4).*exp( - (a - x(1)).^2/(2*sigma_xy^2) - (b - x(2)).^2/(2*sigma_xy^2) - (c - x(3)).^2/(2*sigma_z^2)) + x(5);
R = corrcoef( real_fun(vals), int_array);
    
%Output structure
fit.sigma_y  = sigma_xy;
fit.sigma_x  = sigma_xy;
fit.sigma_z  = sigma_z;
fit.y_pos    = vals(1);
fit.x_pos    = vals(2);
fit.z_pos    = vals(3);
fit.int      = vals(4);
fit.gauss_bg = vals(5);
fit.dist_2_periphery = NaN;
fit.cell_id = NaN;
fit.resnorm = resnorm;
fit.corrcoef = R(1,2);

%Raw signal. Estimate the local background from a MASK defined by guass
%params. 

    %First find distance from center to all other input data points from
    %stack
    D = (( fit.x_pos - b ).^2 + ( fit.y_pos - a ).^2 + ( fit.z_pos - c ).^2 ) .^0.5;
    %Find the entries between params.bg_area(2) and (1)
    sel = D >= params.bg_mask(1) & D <= params.bg_mask(2);
    %Get the mean image intensity of pixels within these bounds
    fit.real_bg = mean( int_array( sel ) );

    %Now find points that are less than the bg_mask(1) boundary
    sel = D <= params.bg_mask(1);
    %Sum up pixels above background
    fit.raw_sum = sum( int_array( sel ) - fit.real_bg );


end




