%% Use imclose
function BW = close_filter(BW, imclose_r)
  
    stats = regionprops(BW,'Centroid','Area','PixelList','PixelIdxList');

    %Now use imclose to remove gaps. Operate on individual cells so
    %that we minimize connecting cells that are close together. 
    se= strel('disk',imclose_r,8);
    new_BW= zeros(size(BW));


    for c = 1:length(stats)
        X = stats(c).PixelList(:,1);
        Y = stats(c).PixelList(:,2);

        x_range = [min(X):max(X)];
        y_range = [min(Y):max(Y)];

        %Need to make a sub-img! 
        sub_img = false(length(y_range),length(x_range));
        %Shift og pixels to sub_img pixels. 
        X = X - x_range(1) + 1;
        Y = Y - y_range(1) + 1;   
        %Get index. 
        ind = sub2ind(size(sub_img),Y,X);
        sub_img(ind) = 1;
        %Close. 
        sub_img = imclose(sub_img,se);

        %Find closed coordinates. 
        [Y,X] = ind2sub(size(sub_img),find(sub_img));
        %Shift back into place. 
        X = X + x_range(1) - 1;
        Y = Y + y_range(1) - 1;

        %New index. 
        ind = sub2ind(size(BW),Y,X);

        %Fill in the large image. 
        new_BW(ind) = 1;            

    end
    BW = logical(new_BW);

end