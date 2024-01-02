% converts a pixelidxlist to a vector 
function p = index_2_vector(img_size, idx)

switch length(img_size)

    case 2
        
        [pY,pX] = ind2sub(img_size, idx);
        p = [pX,pY];

    case 3
    
        [pY,pX,pZ] = ind2sub(img_size, idx);
        p = [pX,pY,pZ];

end