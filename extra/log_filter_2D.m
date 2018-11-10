%Implementing separable 3D log filter

function img_filt = log_filter_2D( img, sigma, filter_size)

%If the filter_size isn't supplied, calculate it:
if(nargin<3)
    filter_size = 2*floor(4.*sigma + 0.5) + 1;
end

%Sigma is vector. Should be sigmax, sigmay, sigmaz. 
sigma_x = sigma(1);
sigma_y = sigma(2);

%Generate the x,y,z values for the LoG filter kernel. 
x = [1:filter_size(1)];
y = [1:filter_size(2)];

%Now center
x = x - mean(x);
y = y - mean(y);


%% LoG filter done as the sum of three separate filters: 

%Integration constant
C = 1 / ( (2*pi)^1.5*sigma_x*sigma_y*sigma_z);
%Base functions
h_f = @(x,sigma) C*( (x.^2./sigma^4) - 1/sigma^2 ).*exp(-x.^2./( 2*sigma^2));
h_b = @(x,sigma) exp(-x.^2/(2*sigma^2));

%First pre-compute the individual detectors
h_f_x = h_f(x,sigma_x);
h_f_y = h_f(y,sigma_y);
h_f_z = h_f(z,sigma_z);

h_b_x = h_b(x,sigma_x);
h_b_y = h_b(y,sigma_y);
h_b_z = h_b(z,sigma_z);

%Do some reshaping
h_f_x = reshape(h_f_x,[1, length(h_f_x), 1]);
h_f_y = reshape(h_f_y,[length(h_f_x), 1, 1]);
h_f_z = reshape(h_f_z,[1, 1, length(h_f_z)]);

h_b_x = reshape(h_b_x,[1, length(h_b_x), 1]);
h_b_y = reshape(h_b_y,[length(h_b_y), 1, 1]);
h_b_z = reshape(h_b_z,[1, 1, length(h_b_z)]);

%% Now convolve the image. 

%Individual separable filters
I_1 = imfilter( imfilter( imfilter( img, h_f_x,'same','replicate'), h_b_y,'same','replicate'),h_b_z,'same','replicate');
I_2 = imfilter( imfilter( imfilter( img, h_b_x,'same','replicate'), h_f_y,'same','replicate'),h_b_z,'same','replicate');
I_3 = imfilter( imfilter( imfilter( img, h_b_x,'same','replicate'), h_b_z,'same','replicate'),h_f_z,'same','replicate');

img_filt = I_1 + I_2 + I_3;

end



