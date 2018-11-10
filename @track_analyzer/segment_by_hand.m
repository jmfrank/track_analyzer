%Edit a segmented image by hand-drawing boundaries. 
%1-25-18: written to quickly make nuclear and cyto ROIs for high-frequency
%imaging where the roi's can be constant for a time series. 
%Expecting 2-channel images 
function obj = segment_by_hand(obj,params)

    %Defaults. 
    params = default_params( params );
    
    
    daspect([1,1,1])

    %clear the field
    obj.masks = [];
    
    debug =0;
        
    %Look too see if a max-p image exists. 
    if( exist( obj.exp_info.max_p_img,'file' ) )


        %Just one image. 
        fname = obj.exp_info.max_p_img;

        [reader,X,Y,Z,C,T] = bfGetReader(fname);
        this_img = bfopen(fname);

        frame_range = [1:T]';
        img_val = 1.*ones(length(frame_range),1);
        %Add to frame2img. 
        frame2img=[frame_range,frame_range, img_val];

        Img = zeros(Y,X,T);
        for t = 1:T
            plane = get_planesZCT(reader,Z,params.view_channel,t);
            Img(:,:,t) = this_img{1}{plane,1};
        end
        %Close reader
        reader.close();
    else
        error('missing max-p img');
        
    end
    
    %First display. 
    figure(1);
    clf
    imshow3D( Img );
    pause %< pause so we can scan to see how many cells were bleached. 
    
    %Ask for the number of cells: 
    prompt = {'Enter the number of cells:'};
    dlgtitle = 'No. cells';
    numlines = 1;
    defaultans = {'1'};
    answer = inputdlg(prompt,dlgtitle,numlines,defaultans);
    
    n_cells = str2num( answer{:} );
    
    %Use a while loop to allow repeating if mistake was made
    pass_seg = 0;
    
    hold on

    while(~pass_seg)

        %Loop over the number of cells
        for c = 1:n_cells

            %First create nuclear ROI mask;
            H = imfreehand('Closed',1);

            nuc_mask = H.createMask;
            %imagesc( nuc_mask )
            
            %Now go back to subplot 1 and create cyto ROI mask
            H = imfreehand('Closed',1);

            cyto_mask = H.createMask;

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
    
    %Auto save. 
    obj.save;
    
    


end
