%% Fill holes
function BW = fill_holes(BW)

    if length(size(BW))==2
        BW = imfill(BW,'holes');
    elseif length(size(BW))==3

        pad_size = 60;
        Z = size(BW,3);
        %Now fill in some holes. 
        for i = 1:Z
            %Also fill holes laterally. Vertical holes can be
            %problematic if we just to imfill with 3D stack

            %first pad 
            pad_plane = padarray(BW(:,:,i),[pad_size,pad_size],'symmetric','both');
            pad_plane_fill = imfill(pad_plane,'holes');

            BW(:,:,i) = pad_plane_fill(pad_size+1:end-pad_size, pad_size+1: end-pad_size);
        end
    
    end
end