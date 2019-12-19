%% diffusion gradients
function J = diff_gradients(I,params)
    
    %%  Step1: Denoising Filters
    %G(:,:,zslice) = imfilter(I(:,:,zslice), h,'replicate');
    %%Optional filters median or deconvlucy
    
    %% Check dimensionality. 
    dims = length(size(I));
    
    switch dims
        
        case 2

            G = medfilt2(I, params.median_filter(1:2));
            %G(:,:,zslice) = deconvlucy(G(:,:,zslice),h); %%use the same 'h' as for Gaussian or define new variable

            %%Step2: Perona & Malik non-linear isotropic diffusion filter, refer to
            %%diffusioncode.m for details, here alter only diffuse_iterations & kappa1       
            Diff_im = diffusioncode(I, params.diffuse_iterations, 0.1429, params.kappa1, params.kappa2, params.option);
   
            %%Step3a: Gauss gradient, 1st derivative of the image, refer to gaussgradient.m for details
            Fim=mat2gray(Diff_im);
            [imx,imy]=gaussgradient(Fim,params.sigmagradient(1));
            Mag_fim = hypot(imx,imy);

            %%Step3b: Laplacian
            [L_imxx, L_imxy] = gaussgradient(imx,params.sigmagradient(1));
            [L_imyx, L_imyy] = gaussgradient(imy,params.sigmagradient(1));

        case 3
            G = zeros(size(I));
            Diff_im = zeros(size(I));
            %Median filtering and diffusion in 2D. 
            for i  = 1:size(I,3)
                G(:,:,i) = medfilt2(I(:,:,i),params.median_filter(1:2));
                Diff_im(:,:,i) = diffusioncode(I(:,:,i), params.diffuse_iterations, 0.1429, params.kappa1, params.kappa2, params.option);
            end
            
            %Rescale to [0,1]
            Fim = mat2gray(Diff_im); clear Diff_im;
            %Get gauss gradient. Using separable filters. sigmagradient now 3 component vector
            [imx,imy]=gaussgradient2D_sep(Fim,params.sigmagradient); clear Fim;
            %Total magnitude. 
            Mag_fim = hypot(imx,imy);

            %Laplacian. 
            [L_imxx, L_imxy] = gaussgradient2D_sep(imx,params.sigmagradient);
            [L_imyx, L_imyy] = gaussgradient2D_sep(imy,params.sigmagradient);

    end

    Mag_Lim = L_imxx+L_imyy;  %Laplacian (trace of the 2D matrix)
    Mag_Lim=Mag_Lim.*(Mag_Lim>0);

    %%Step3c: Hessian Determniant
    Det_hessian=L_imxx.*L_imyy-L_imxy.*L_imyx;
    Det_hessian=Det_hessian.*(Det_hessian<0);
    %%Step4: Masking function using tanh on the Summed Derivatives from Step3
    X=params.gamma-((params.alpha*Mag_fim + params.beta*Mag_Lim + params.epsilon*abs(Det_hessian)) /params.delta);
    Multi=0.5*(tanh(X)+1);% masking function
    J=(double(G).*Multi); % masked image applied on smoothened image
    
end