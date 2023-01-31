%Fit a 3D gaussian to 3D intensity data

%Input 3D image stack
%Centroid is [y,x] estimated center of blob
% y,x,z are the 


function fit = fit_2D_gauss( img, centroid, params, FG)

%x is the fit to the 3D gaussian parameters are: 
%    x(1) = sigma-x/y
%    x(2) = sigma-z
%    x(3) = y-location
%    x(4) = x-location
%    x(5) = z-location
%    x(6) = scaling factor
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
fun = @(x) x(4).*exp( - (a - x(2)).^2./(2*x(1)^2) - (b - x(3)).^2./(2*x(1)^2)) + x(5) - int_array;

%Perform non-linear least squares:

    %Set up initial guess:
    x0(1) = params.xy_sigma_est;
    x0(2) = centroid(2);
    x0(3) = centroid(1);
    x0(4) = params.int_guess;
    x0(5) = params.bg_guess;
    
    lb = [0.2,1,1,1,0];
    ub = [10,size(img,1),size(img,2),5000,500];
    options = optimoptions('lsqnonlin','Display','off');
    options.Algorithm = 'levenberg-marquardt';

    [vals, resnorm] = lsqnonlin(fun,x0,[],[],options);
 
%Correlation-coefficient
real_fun = @(x) x(4).*exp( - (a - x(2)).^2./(2*x(1)^2) - (b - x(3)).^2./(2*x(1)^2)) + x(5)
R = corrcoef( real_fun(vals), int_array);
    
%Output structure
fit.sigma       = [vals(1),vals(1)]; %<remember y,x,z
fit.pos         = vals(2:3); %<remember y,x,z
fit.int         = vals(4);
fit.gauss_bg    = vals(5);
fit.dist_2_periphery = NaN;
fit.cell_id = NaN;
fit.resnorm = resnorm;
fit.corrcoef = R(1,2);

%Raw signal. 
    %Find the entries between params.bg_area(2) and (1)
    %Get the mean image intensity of pixels within these bounds
    fit.real_bg =params.bg_guess;
    %Now find points that are less than the bg_mask(1) boundary
    %Sum up pixels above background
    fit.raw_sum = sum( int_array(FG(:)) - fit.real_bg );
    
% figure(7);
% imagesc(img)
% 
% hold on
% plot(fit(3),fit(2),'*r')
% hold off

end




