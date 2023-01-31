%Fit a 3D gaussian to 3D intensity data

%Input 3D image stack
%Centroid is [y,x] estimated center of blob
% y,x,z are the 

function [fit,out_stack] = fit_3D_gauss_SIGxy_constant( stack, centroid, params, FG)

%x is the fit to the 3D gaussian parameters are: 
%    x(1) = sigma-xy
%    x(2) = sigma-z
%    x(3) = y-location
%    x(4) = x-location
%    x(5) = z-location
%    x(6) = scaling factor
%    x(7) = background
%    a = y 
%    b = x 
%    c = z 

%Linearize stack
int_array = stack(:);
%Get an array list for the coordinates
[I,J,K] = ind2sub(size(stack),[1:length(int_array)]);
a = I';
b = J';
c = K';

%Create 3D gaussian function to minimize:
fun = @(x) x(6).*exp( - (a - x(3)).^2/(2*x(1)^2) - (b - x(4)).^2/(2*x(1)^2) - (c - x(5)).^2/(2*x(2)^2)) + x(7) - int_array;

%Perform non-linear least squares:

    %Set up initial guess:
    x0(1) = params.xy_sigma_est; %Guess y and x sigma to be the same initially
    x0(2) = params.z_sigma_est;
    x0(3) = centroid(2);
    x0(4) = centroid(1);
    x0(5) = centroid(3); %Guess halfway into stack
    x0(6) = params.int_guess;
    x0(7) = params.bg_guess;
    
    %Lower bound
    pos_lb = size(stack)*0.25;
    pos_ub = size(stack)*0.75;
    lb = [0.5,0.5,pos_lb,1,0.2*params.bg_guess];
    %Upper bound. Don't let location move very far from centroid. 
    ub = [2,  2,  pos_ub,65000,1.5*params.bg_guess];
    
    
    options = optimoptions('lsqnonlin','Display','off');
    %options.Algorithm = 'levenberg-marquardt';
    [vals, resnorm] = lsqnonlin(fun,x0,lb,ub,options);
    
%Correlation-coefficient
real_fun = @(x) x(6).*exp( - (a - x(3)).^2/(2*x(1)^2) - (b - x(4)).^2/(2*x(1)^2) - (c - x(5)).^2/(2*x(2)^2)) + x(7);
R = corrcoef( real_fun(vals), int_array);

%Output structure
fit.sigma       = [vals(1),vals(1),vals(2)]; %<remember y,x,z
fit.pos         = vals(3:5); %<remember y,x,z
fit.int         = vals(6);
fit.gauss_bg    = vals(7);
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
    
%Output a theoretical image. 
if(nargout==2)
   out_int = real_fun(vals);
   out_stack = reshape(out_int,size(stack));
    
end

end




