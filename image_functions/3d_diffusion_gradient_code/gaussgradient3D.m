function [gx,gy,gz]=gaussgradient3D(IM,sigma)
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

%generate a 3-D Gaussian kernel along x direction
hx = zeros(size);
for i=1:size(1)
    for j=1:size(2)
        for k = 1:size(3)
            %Coordinate. 
            u=[i-halfsize-1 j-halfsize-1 k-halfsize-1];
            %Gauss gradient. 
            hx(i,j,k)=gauss(u(1),sigma(1))*gauss(u(3),sigma(3))*dgauss(u(2),sigma(2));
        end
    end
end
%Normalize. 
sq_hx = abs(hx).^2;
norm_val = sqrt( sum( sq_hx(:)));
hx = hx ./ norm_val;
%hx=hx/sqrt(sum(sum(abs(hx).*abs(hx))));

%generate a 3-D Gaussian kernel along y direction
if(sigma(1)==sigma(2))
    hy = permute(hx,[2,1,3]);
else
    error('missing some code here')
end

%Now we need to generate 3D gaussian kernel in z direction. 
hz = zeros(size);
for i=1:size(1)
    for j=1:size(2)
        for k = 1:size(3)
            %Coordinate. 
            u=[i-halfsize-1 j-halfsize-1 k-halfsize-1];
            %Gauss gradient. 
            hz(i,j,k)=gauss(u(1),sigma(1))*gauss(u(2),sigma(2))*dgauss(u(3),sigma(3));
        end
    end
end

%3-D filtering
gx=imfilter(IM,hx,'replicate','conv');
gy=imfilter(IM,hy,'replicate','conv');
gz=imfilter(IM,hz,'replicate','conv');

function y = gauss(x,sigma)
%Gaussian
y = exp(-x^2/(2*sigma^2)) / (sigma*sqrt(2*pi));
function y = dgauss(x,sigma)
%first order derivative of Gaussian
y = -x * gauss(x,sigma) / sigma^2;