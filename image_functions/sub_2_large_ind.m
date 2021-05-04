% Convers the index of matrices from a sub-region to the original matrix.
% Important for performing tasks on smaller part of image.

function idx = sub_2_large_ind( sub_ind, ranges, sub_size, og_size )
    
    switch length( og_size )

        case 2

            [y,x] = ind2sub(sub_size,sub_ind);

            % Shift to large image. 
            X = x + ranges.x(1) -1;
            Y = y + ranges.y(1) -1;

            % Now get index within large image. 
            idx = sub2ind(og_size, Y,X);
            
        case 3
            
            [y,x,z] = ind2sub(sub_size,sub_ind);

            % Shift to large image. 
            X = x + ranges.x(1) -1;
            Y = y + ranges.y(1) -1;
            Z = z + ranges.z(1) -1;

            % Now get index within large image. 
            idx = sub2ind(og_size, Y,X,Z);

    end
end