function [gx,gy,gz]=gaussgradient2D_sep(IM,sigma)
%GAUSSGRADIENT Gradient using first order derivative of Gaussian.
%  [gx,gy, gz]=gaussgradient(IM,sigma) outputs the gradient image gx, gy, gz of
%  image IM using a 3-D Gaussian kernel. Sigma is the standard deviation of
%  this kernel along each direction.
%
%  Contributed by Guanglei Xiong (xgl99@mails.tsinghua.edu.cn)
%  at Tsinghua University, Beijing, China.

%determine the appropriate size of kernel. The smaller epsilon, the larger
%size.
epsilon=1e-2;
halfsize=ceil(sigma.*sqrt(-2.*log(sqrt(2*pi).*sigma.*epsilon)));
size=2*halfsize+1;

% First generate a 3-D Gaussian kernel along x direction
hx = zeros(size(1:2));
for i=1:size(1)
    for j=1:size(2)
        %Coordinate. 
        u=[i-halfsize(1)-1 j-halfsize(2)-1];
        %Gauss gradient. 
        hx(i,j)=gauss(u(1),sigma(1))*dgauss(u(2),sigma(2));
    end
end
%Normalize. 
sq_hx = abs(hx).^2;
norm_val = sqrt( sum( sq_hx(:)));
hx = hx ./ norm_val;
%Decompose. 
[K,~,err] = SeparateKernel(hx);
%Filter in each direction. 
try
    gx = gather(imfilter(imfilter(gpuArray(IM),K{1},'replicate','conv'),K{2},'replicate','conv'));
catch
    gx = imfilter(imfilter(IM,K{1},'replicate','conv'),K{2},'replicate','conv');
end

%Now switch filter to y direction. generate a 3-D Gaussian kernel along y direction
if(sigma(1)==sigma(2))
    hy = hx';
else
    error('missing some code here')
end
[K,~,err] = SeparateKernel(hy);
try
    gy = gather(imfilter(imfilter(gpuArray(IM),K{1},'replicate','conv'),K{2},'replicate','conv'));
catch
    gy = imfilter(imfilter(IM,K{1},'replicate','conv'),K{2},'replicate','conv');
end




function y = gauss(x,sigma)
%Gaussian
y = exp(-x^2/(2*sigma^2)) / (sigma*sqrt(2*pi));
function y = dgauss(x,sigma)
%first order derivative of Gaussian
y = -x * gauss(x,sigma) / sigma^2;