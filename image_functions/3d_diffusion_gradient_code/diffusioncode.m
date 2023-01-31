% Adapated by Bhavna Rajasekaran

% References:
%   P. Perona and J. Malik.
%   Scale-Space and Edge Detection Using Anisotropic Diffusion.
%   IEEE Transactions on Pattern Analysis and Machine Intelligence,
%   12(7):629-639, July 1990.
%   P. D. Kovesi. MATLAB and Octave Functions for Computer Vision and Image Processing.
%   School of Computer Science & Software Engineering,
%   The University of Western Australia. Available from:
%   <http://www.csse.uwa.edu.au/~pk/research/matlabfns/>.
% 
% Original Code: anisodiff3D.m
% Daniel Simoes Lopes
% ICIST
% Instituto Superior Tecnico - Universidade Tecnica de Lisboa
% danlopes (at) civil ist utl pt
% http://www.civil.ist.utl.pt/~danlopes

% Code modified to include forward and backward anisotropic diffusion
%On the combined forward and backward anisotropic diffusion
%scheme for the multispectral image enhancement
%Bogdan Smolka Rastislav Lukac

function diff_im = diffusioncode(im, num_iter, delta_t, kappa1,kappa2, option)
%ANISODIFF2D Conventional anisotropic diffusion
%   DIFF_IM = ANISODIFF2D(IM, NUM_ITER, DELTA_T, KAPPA, OPTION) perfoms
%   conventional anisotropic diffusion (Perona & Malik) upon a gray scale
%   image. A 2D network structure of 8 neighboring nodes is considered for
%   diffusion conduction.
%
%       ARGUMENT DESCRIPTION:
%               IM       - gray scale image (MxN).
%               NUM_ITER - number of iterations.
%               DELTA_T  - integration constant (0 <= delta_t <= 1/7).
%                          Usually, due to numerical stability this
%                          parameter is set to its maximum value.
%               KAPPA    - gradient modulus threshold that controls the conduction.
%               OPTION   - conduction coefficient functions proposed by Perona & Malik:
%                          1 - c(x,y,t) = exp(-(nablaI/kappa).^2),
%                              privileges high-contrast edges over low-contrast ones.
%                          2 - c(x,y,t) = 1./(1 + (nablaI/kappa).^2),
%                              privileges wide regions over smaller ones.

%%% Extra options for forward-backward diffusions%%%%%
%                          3- Forward-backward diffusion c(x,y,t) =
%                          2.*exp(-(nablaI/kappa1).^2)-exp(-(nablaI/kappa2).^2),
%                          4- Forward-backward diffusion c(x,y,t) =2./(1 +(nablaI/kappa1).^2)-1./(1 + (nablaI/kappa2).^2)
%       OUTPUT DESCRIPTION:
%                DIFF_IM - (diffused) image with the largest scale-space parameter.
%
%

% Convert input image to double.
im = double(im);

% PDE (partial differential equation) initial condition.
diff_im = im;

% Center pixel distances.
dx = 1;
dy = 1;
%dd = sqrt(2);
%dx = 0.5;
%dy = 0.5;
dd = sqrt(dx^2+dy^2);

% 2D convolution masks - finite differences.
hN = [0 1 0; 0 -1 0; 0 0 0];
hS = [0 0 0; 0 -1 0; 0 1 0];
hE = [0 0 0; 0 -1 1; 0 0 0];
hW = [0 0 0; 1 -1 0; 0 0 0];
hNE = [0 0 1; 0 -1 0; 0 0 0];
hSE = [0 0 0; 0 -1 0; 0 0 1];
hSW = [0 0 0; 0 -1 0; 1 0 0];
hNW = [1 0 0; 0 -1 0; 0 0 0];

% Anisotropic diffusion.
    for t = 1:num_iter
    
    % Finite differences. [imfilter(.,.,'conv') can be replaced by conv2(.,.,'same')]
    nablaN = imfilter(diff_im,hN,'conv');
    nablaS = imfilter(diff_im,hS,'conv');
    nablaW = imfilter(diff_im,hW,'conv');
    nablaE = imfilter(diff_im,hE,'conv');
    nablaNE = imfilter(diff_im,hNE,'conv');
    nablaSE = imfilter(diff_im,hSE,'conv');
    nablaSW = imfilter(diff_im,hSW,'conv');
    nablaNW = imfilter(diff_im,hNW,'conv');
    
    % Diffusion function.
    if option == 1
        kappa2=0;
        cN = exp(-(nablaN/kappa1).^2)-kappa2;
        cS = exp(-(nablaS/kappa1).^2)-kappa2;
        cW = exp(-(nablaW/kappa1).^2)-kappa2;
        cE = exp(-(nablaE/kappa1).^2)-kappa2;
        cNE = exp(-(nablaNE/kappa1).^2)-kappa2;
        cSE = exp(-(nablaSE/kappa1).^2)-kappa2;
        cSW = exp(-(nablaSW/kappa1).^2)-kappa2;
        cNW = exp(-(nablaNW/kappa1).^2)-kappa2;
    elseif option == 2
        kappa2=0;
        cN = 1./(1 + (nablaN/kappa1).^2)-kappa2;
        cS = 1./(1 + (nablaS/kappa1).^2)-kappa2;
        cW = 1./(1 + (nablaW/kappa1).^2)-kappa2;
        cE = 1./(1 + (nablaE/kappa1).^2)-kappa2;
        cNE = 1./(1 + (nablaNE/kappa1).^2)-kappa2;
        cSE = 1./(1 + (nablaSE/kappa1).^2)-kappa2;
        cSW = 1./(1 + (nablaSW/kappa1).^2)-kappa2;
        cNW = 1./(1 + (nablaNW/kappa1).^2)-kappa2;
    elseif option == 3
        cN = 2.*(exp(-(nablaN/kappa1).^2))-exp(-(nablaN/kappa2).^2);
        cS = 2.*(exp(-(nablaS/kappa1).^2))-exp(-(nablaS/kappa2).^2);
        cW = 2.*(exp(-(nablaW/kappa1).^2))-exp(-(nablaW/kappa2).^2);
        cE = 2.*(exp(-(nablaE/kappa1).^2))-exp(-(nablaE/kappa2).^2);
        cNE = 2.*(exp(-(nablaNE/kappa1).^2))-exp(-(nablaNE/kappa2).^2);
        cSE = 2.*(exp(-(nablaSE/kappa1).^2))-exp(-(nablaSE/kappa2).^2);
        cSW = 2.*(exp(-(nablaSW/kappa1).^2))-exp(-(nablaSW/kappa2).^2);
        cNW = 2.*(exp(-(nablaNW/kappa1).^2))-exp(-(nablaNW/kappa2).^2);
    elseif option == 4
        cN = 2./(1 + (nablaN/kappa1).^2)-1./(1 + (nablaN/kappa2).^2);
        cS = 2./(1 + (nablaS/kappa1).^2)-1./(1 + (nablaS/kappa2).^2);
        cW = 2./(1 + (nablaW/kappa1).^2)-1./(1 + (nablaW/kappa2).^2);
        cE = 2./(1 + (nablaE/kappa1).^2)-1./(1 + (nablaE/kappa2).^2);
        cNE = 2./(1 + (nablaNE/kappa1).^2)-1./(1 + (nablaNE/kappa2).^2);
        cSE = 2./(1 + (nablaSE/kappa1).^2)-1./(1 + (nablaSE/kappa2).^2);
        cSW = 2./(1 + (nablaSW/kappa1).^2)-1./(1 + (nablaSW/kappa2).^2);
        cNW = 2./(1 + (nablaNW/kappa1).^2)-1./(1 + (nablaNW/kappa2).^2);
    else
        error('Incorrect option');
    end
    
    % Discrete PDE solution.
    diff_im = diff_im + ...
        delta_t*(...
        (1/(dy^2))*cN.*nablaN + (1/(dy^2))*cS.*nablaS + ...
        (1/(dx^2))*cW.*nablaW + (1/(dx^2))*cE.*nablaE + ...
        (1/(dd^2))*cNE.*nablaNE + (1/(dd^2))*cSE.*nablaSE + ...
        (1/(dd^2))*cSW.*nablaSW + (1/(dd^2))*cNW.*nablaNW );
 
    end
end