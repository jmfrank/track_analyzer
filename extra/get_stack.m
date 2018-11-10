%Quick way to get a z-stack of a bf image reader at time t and channel c. 
function stack = get_stack( reader, t, c)


%Get the image size of this series. 
size_x = reader.getSizeX;
size_y = reader.getSizeY;
Z = reader.getSizeZ;

%Get the images for this z-stack according to bioformats reader
stack = zeros(size_y,size_x,Z);
for i = 1:Z
    this_plane = reader.getIndex(i-1,c-1,t-1)+1;
    stack(:,:,i) = bfGetPlane(reader,this_plane);
end


end
