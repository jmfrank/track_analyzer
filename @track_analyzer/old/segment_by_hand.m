%Edit a segmented image by hand-drawing boundaries. 
%1-25-18: written to quickly make nuclear and cyto ROIs for high-frequency
%imaging where the roi's can be constant for a time series. 
%Expecting 2-channel images 
function obj = segment_by_hand(obj)

    daspect([1,1,1])

    %clear the field
    obj.masks = [];
    
    debug =0;
    
    Z = obj.exp_info.z_frames;
    T = obj.exp_info.t_frames;

    %First plot frame 1. 
        t=1;
        start = Z.*(t-1)*2 + 1;
        indices  = start:2:start + 2*Z-1;
        %Grab this Z-stack
        IMG = tiffread( obj.exp_info.img_file, indices );

        %Detect number of channels
        n_channels = IMG(1).lsm.DimensionChannels;
        
        if(n_channels==1)

            %Getting size of example image
            [y,x] = size( IMG(1).data );

            %Now collect data for each channel. Just use channel 1 for both images. 
            img{1} = zeros(y,x,Z);
            img{2} = zeros(y,x,Z);
            for z = 1:Z
              img{1}(:,:,z) = IMG(z).data;
              img{2}(:,:,z) = IMG(z).data;
            end          
            
        elseif( n_channels == 2)

            %Getting size of example image
            [y,x] = size( IMG(1).data{1} );

            %Now collect data for each channel
            img{1} = zeros(y,x,Z);
            img{2} = zeros(y,x,Z); 
            for z = 1:Z
              img{1}(:,:,z) = IMG(z).data{1};
              img{2}(:,:,z) = IMG(z).data{2};
            end
        end
        
        mip{1} = max(img{1},[],3);

        subplot(2,2,1)
        imagesc(mip{1})
        colormap('gray');

    %Now plot the last frame
        t=T;
        start = Z.*(t-1)*2 + 1;
        indices  = start:2:start + 2*Z-1;
        %Grab this Z-stack
        IMG = tiffread( obj.exp_info.img_file, indices );

        if(n_channels==1)

            %Getting size of example image
            [y,x] = size( IMG(1).data );

            %Now collect data for each channel. Just use channel 1 for both images. 
            img{1} = zeros(y,x,Z);
            img{2} = zeros(y,x,Z);
            for z = 1:Z
              img{1}(:,:,z) = IMG(z).data;
              img{2}(:,:,z) = IMG(z).data;
            end          
            
        elseif( n_channels == 2)

            %Getting size of example image
            [y,x] = size( IMG(1).data{1} );

            %Now collect data for each channel
            img{1} = zeros(y,x,Z);
            img{2} = zeros(y,x,Z); 
            for z = 1:Z
              img{1}(:,:,z) = IMG(z).data{1};
              img{2}(:,:,z) = IMG(z).data{2};
            end
        end
       

        mip{1} = max(img{1},[],3);
        mip{2} = max(img{2},[],3);
        
        subplot(2,2,2)
        imagesc(mip{1})
        colormap('gray')
        subplot(2,2,3)
        imagesc(mip{2})
        colormap('gray')
    
    
    %Ask for the number of cells: 
    prompt = {'Enter the number of cells:'};
    dlgtitle = 'No. cells';
    numlines = 1;
    defaultans = {'1'};
    answer = inputdlg(prompt,dlgtitle,numlines,defaultans);
    
    n_cells = str2num( answer{:} );
    
    %Use a while loop to allow repeating if mistake was made
    pass_seg = 0;
    
    while(~pass_seg)

        %Loop over the number of cells
        for c = 1:n_cells

            subplot(2,2,1);

            %First create nuclear ROI mask;
            H = imfreehand('Closed',1);

            nuc_mask = H.createMask;
            subplot(2,2,4);
            hold on
            imagesc( nuc_mask )
            hold off
            
            %Now go back to subplot 1 and create cyto ROI mask
            subplot(2,2,1)
            H = imfreehand('Closed',1);

            cyto_mask = H.createMask;
            subplot(2,2,4);
            hold on
            imagesc( cyto_mask ) 
            hold off

            %Add the masks to the obj
            obj.masks(c).nuc_mask = nuc_mask;
            obj.masks(c).cyto_mask = cyto_mask;

        end
        
        %Use prompt to ask if things went well
        button = questdlg('Pass to move forward, OR re-do','Question','Pass','REDO','Pass');
        
        switch button
            case 'Pass'
                pass_seg = 1;
            case 'REDO'
                pass_seg = 0;
        end
    end
    
    
    


end
