% Generate list of [Y,X,Z] coordinates of a mask. 

function [Y,X,Z] = get_mask_coordinates( mask )


    % Check dimensionality. 
    dims = length(size(mask));
    
    switch dims
        
        
        case 2
            
            [Y,X] = find( logical(mask) );
            
            Z = [];
        case 3
            
            idx = find( logical(mask) );
            
            [Y,X,Z] = ind2sub( size(mask), idx );
            
    end
    



end