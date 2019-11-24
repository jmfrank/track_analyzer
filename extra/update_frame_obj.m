% Simple script to update frame_obj to new format. Will assume channel 1
% segmentation if one channel. Or channel 2 segmentation if 2 channels. 

function update_frame_obj( info )


        obj = get_exp(info);
        
        
        
        if obj.exp_info.channels == 1
            
            seg_channel=1;
            
        elseif obj.exp_info.channels == 2
            
            warning('multi-channel image. assuming channel 2 segmented')
            seg_channel = 2;
            
        end
        
        % Define channel string. 
        channel_str = ['seg_channel_',pad(num2str(seg_channel),2,'left','0')];

        F = obj.get_frame_files;
        
        %loop over frames. 
        for i = 1:length(F)
            
            
            D = load(F{i},'frame_obj');
            
            % Skip if updated. 
            if isfield(D.frame_obj,channel_str)
                warning([F{i}, ' already updated'])
                continue
            end
            
            new_frame_obj = rmfield(D.frame_obj,{'PixelIdxList','centroids','contours','BW'});

            new_struct = struct();
            new_struct.PixelIdxList = D.frame_obj.PixelIdxList;
            new_struct.centroids = D.frame_obj.centroids;
            new_struct.contours = D.frame_obj.contours;
            new_struct.BW = D.frame_obj.BW;
            
            new_frame_obj.(channel_str) = new_struct;
            
            %Now save. 
            parsave(F{i},new_frame_obj)
            
            
        end
end


%Parallel saving technique
function parsave(fname, frame_obj)
try
    save(fname,'frame_obj');
catch
    %For large files...
    save(fname, 'frame_obj','-v7.3');
end
end