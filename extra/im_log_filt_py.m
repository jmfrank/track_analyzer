
%LoG filter implemented to match python code in datadoctor
function img = im_log_filt_py(im, sigma )

    %Calculate a width of the gaussian based on sigma
    w = 2*floor(4*sigma + 0.5) + 1;
    
    %Now generate LoG filter
    h = fspecial('log',w,sigma);
    
    %Now filter image
    img = imfilter(im,h,'symmetric');
    
    
    
end