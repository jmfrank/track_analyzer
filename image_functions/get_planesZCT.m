%% Custom get_planes function. 
function planes = get_planesZCT(reader,Z,C,t)

    for z = 1:Z
        planes(z) = reader.getIndex(z-1,C-1,t-1)+1;
    end


end
