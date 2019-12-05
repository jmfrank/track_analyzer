%Implementing separable 3D log filter

function img_filt = log_filter_2D( img, sigma, filter_size)

%If the filter_size isn't supplied, calculate it:
if(nargin<3)
    filter_size = 2*floor(4.*sigma + 0.5) + 1;
end

%w=fspecial('log', sigma(1:2), filter_size);
%img_filt=imfilter(img,w,'replicate');




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
C = 1 / ( (2*pi)^0.5*sigma_x^2);
%Base functions
h_1 = @(x,sigma) C*( 1 - x.^2/sigma^2).*exp(-x.^2./( 2*sigma^2));
h_2 = @(x,sigma) -C*exp(-x.^2/(2*sigma^2));

%First pre-compute the individual detectors
h1 = h_1(x,sigma_x);
h2 = h_2(x,sigma_x);

%Do some reshaping
h1_x = reshape(h1,[1, length(h1)]);
h1_y = reshape(h1,[length(h1), 1]);

h2_x = reshape(h2,[1,length(h2)]);
h2_y = reshape(h2,[length(h2),1]);

%% Now convolve the image. 
try % utilize gpu. 
    img = gpuArray(img);
    I_1 = imfilter( imfilter(img,h1_x,'same','replicate'),h2_y,'same','replicate');
    I_2 = imfilter( imfilter(img,h2_x,'same','replicate'),h1_y,'same','replicate');
    img_filt = gather(I_1 + I_2);
catch
    I_1 = imfilter( imfilter(img,h1_x,'same','replicate'),h2_y,'same','replicate');
    I_2 = imfilter( imfilter(img,h2_x,'same','replicate'),h1_y,'same','replicate');
    img_filt = I_1 + I_2;
end

end



