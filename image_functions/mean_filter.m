%% Mean Filter
function I = mean_filter(I,MeanFilterSensitivity,MeanFilterNeighborhood,image_bits)

    if image_bits==8
        if ~strcmp(class(I),'uint8')
            I = uint8(I);
        end
        tmp = adaptthresh(I,MeanFilterSensitivity,'NeighborhoodSize',MeanFilterNeighborhood);
    elseif image_bits==12

        %Rescale to 16 bit? 
        I = uint16(I* (2^16-1)/(2^12-1));
        tmp = adaptthresh(I,MeanFilterSensitivity,'NeighborhoodSize',MeanFilterNeighborhood);

    elseif image_bits==16
        if ~strcmp(class(I),'uint16')
            I = uint16(I);
        end
        tmp = adaptthresh(I,MeanFilterSensitivity,'NeighborhoodSize',MeanFilterNeighborhood);

    end

    %%%% This step might help bring everything to 16bit levels. Try
    %%%% manually changing bit size to 16 here. 
    I = tmp.*65000;
end