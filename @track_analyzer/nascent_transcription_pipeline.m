%Pipeline wrapper. Now we open each time point only once, performing
%segmentation / fitting / etc on the frame before continuing. 
function nascent_transcription_pipeline( obj, params, step)

%% Pre-processing using bio-formats. 

%Generate reader. FOR NOW, assume we are looking in series 1. 
reader = bfGetReader(obj.exp_info.img_file);
series = 1;

%Get the number of time points. 
T = reader.getSizeT;

 %Get step.debug
if(step.debug)
    params
    %T = 1;
end  

%% Generate im_info structure
    ZSlicesinStack = obj.exp_info.z_planes;

    image_bits     = reader.getBitsPerPixel;   
    
    %Since we do parallel processing for time points, pre-compute the list
    %of plane indices for each time point. 
    Z_cell = {};
    for t = 1:T
        planes = [];
        for Z = 1:ZSlicesinStack
            planes(Z) = reader.getIndex(Z-1,params.seg_channel-1,t-1)+1;
        end
        Z_cell{t} = planes;
    end
    
%% Go through time points, open the stack, send to various analysis functions as directed by step. 



        
    