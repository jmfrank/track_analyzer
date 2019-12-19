%% CLAHE
function I = clahe(I, tile_size, ClipLimit)     

    % error if 3D. 
    if length(size(I))==3
        
        error('CLAHE function not support in 3D');
        
    end
    
    I_sm = zeros(size(I));
    mu = mean(I(:));
    i_std = std(I(:));
    norm_min = min(I(:));
    norm_max = max(I(:));
    %Set the max val to 4std above mean?
    %norm_max = 4*i_std + mu;
    [size_y,size_x,size_z] = size(I);
    norm_plane = ( I - norm_min ) ./ (norm_max - norm_min );
    tiles = round([size_y /tile_size, size_x / tile_size]);
    I_sm = adapthisteq(norm_plane,'NumTiles',tiles,'ClipLimit',ClipLimit); %,'Distribution','exponential','Alpha',0.1);

    %Replace background 
    %if(step.subtract_bg)
    %    I_sm(sel) = 0;
    %end
    I = I_sm.*65535;
end

% %% 3D functionality?
%     if(step.CLAHE)
%         I_sm = zeros(size(I));
%         %Weights. 
%         lower_bound = min(plane_mean)-0.5.*min(plane_std);
%         max_bound   = max(plane_mean);
%         weights = (plane_mean-lower_bound) ./ (max_bound-lower_bound);
%         for i = 1:ZSlicesinStack
%             this_plane = I(:,:,i);
%             norm_max = max(this_plane(:));
%             norm_min = min(this_plane(:));
%             
%             norm_plane = ( this_plane - norm_min ) ./ (norm_max - norm_min );
%             tiles = round([size_y / params.tile_size, size_x / params.tile_size]);
%             I_sm(:,:,i) = adapthisteq(norm_plane,'NumTiles',tiles,'ClipLimit',params.ClipLimit).*weights(i); %,'Distribution','exponential','Alpha',0.1);
%             
%         end
%         
%         I = I_sm.*65535;
%         clear I_sm
%     end