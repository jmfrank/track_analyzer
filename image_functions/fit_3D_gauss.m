%Fit a 3D gaussian to 3D intensity data

%Input 3D image stack
%Centroid is [y,x] estimated center of blob
% y,x,z are the 

function fit = fit_3D_gauss( stack, centroid, params, FG)

%x is the fit to the 3D gaussian parameters are: 
%    x(1) = sigma-y
%    x(2) = sigma-x
%    x(3) = sigma-z
%    x(4) = y-location
%    x(5) = x-location
%    x(6) = z-location
%    x(7) = scaling factor
%    x(8) = background
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
fun = @(x) x(7).*exp( - (a - x(4)).^2/(2*x(1)^2) - (b - x(5)).^2/(2*x(2)^2) - (c - x(6)).^2/(2*x(3)^2)) + x(8) - int_array;

%Perform non-linear least squares:

    %Set up initial guess:
    x0(1) = params.xy_sigma_est; %Guess y and x sigma to be the same initially
    x0(2) = params.xy_sigma_est; %Guess y and x sigma to be the same initially
    x0(3) = params.z_sigma_est;
    x0(4) = centroid(2);
    x0(5) = centroid(1);
    x0(6) = centroid(3); %Guess halfway into stack
    x0(7) = params.int_guess;
    x0(8) = params.bg_guess;
    
    %Lower bound
    pos_lb = size(stack)*0.25;
    pos_ub = size(stack)*0.75;
    
    lb = [0.5,0.5,0.5,pos_lb,1,0];
    %Upper bound. Don't let location move very far from centroid. 
    ub = [5,  5,  5,  pos_ub,65000,60000];
    options = optimoptions('lsqnonlin','Display','off');
    [vals, resnorm] = lsqnonlin(fun,x0,lb,ub,options);
    
    
%Output structure
fit.sigma       = vals(1:3); %<remember y,x,z
fit.pos         = vals(4:6); %<remember y,x,z
fit.int         = vals(7);
fit.gauss_bg    = vals(8);
fit.dist_2_periphery = NaN;
fit.cell_id = NaN;
fit.resnorm = resnorm;

%Raw signal. Estimate the local background from a MASK defined by guass
%params. 

    %Find the entries between params.bg_area(2) and (1)
    %Get the mean image intensity of pixels within these bounds
    fit.real_bg =params.bg_guess;
    %Now find points that are less than the bg_mask(1) boundary
    %Sum up pixels above background
    fit.raw_sum = sum( int_array(FG(:)) - fit.real_bg );
end




