%% Filter out objects with too big of radius of gyration
function out_stats = filter_rg(stats, max_rg,img_size)
    out_stats=[];
    c=0;
    for i =1:length(stats)

        %Check radius of gyration.
        [y,x] = ind2sub(img_size,stats(i).PixelIdxList);
        %Calculate radius of gryation
        pos=[x,y];
        com=[mean(x),mean(y)];

        rg = sqrt(1/length(y).*sum(sum( (pos-com).^2,2)));

        if rg > max_rg
            continue
        end
        c=c+1;
        out_stats(c)=stats(i);
    end
end