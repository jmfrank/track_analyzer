classdef track_analyzer

    
    properties (SetAccess = private)
        
        
        % The trajectories stored in a cell array, one T x n_dim per particle
        tracks = {};
        % Trajectories corresponding to the spots localized in cell tracks
        spot_tracks = {};
        
        %Cells (originally, 2D contour traces). Need to change for 3D
        %segmentation. 
        cells
        
        %Background value. Each experiment should have a background value
        %in order to calculate intensity properly. 
        background 
        
        %Frap field in case we need to use a frap object. 
        frap
        
        %Signal correlation. 
        sigcorr
        
        %Masks for simple roi stuff
        masks
    end
    
    properties (SetAccess = public )
        results = []; %Results from nascent RNA spot localization
        nuc_cyto_data = [];
        tags = []; %Tags for identifying bleached tracks, bad tracks, etc. 
        
        %Store the experimental info (i.e. file directories, experimental parameters, etc. 
        exp_info = {};
        
    end
     
    
    %% Constructor
    methods
        
        %Simply makes the new obj
        function obj = track_analyzer(exp_info)
            obj.exp_info = exp_info;
        end
        
    end
        
    
    %% Static and public functions
    
    methods
        
        %Save
        function obj = save(obj)

            %Save myself! Use append in case there are other parts to the track_file. 
            track_obj = obj;
            save(obj.exp_info.track_file,'track_obj');

        end


        %Can also load data into the object
        function obj = add_data(obj,in_data)
            
            %Add new tracks to the object
            obj.tracks = in_data.tracks;
            obj.cells  = in_data.cells;
            obj.img_files = in_data.img_files;
        end
        
        function obj = update_field(obj, fieldname, data )
            
            obj.(fieldname) = data;
            
        end
        
        %Clear out cell / track flags. 
        function obj = clear_flags(obj)
            
            try 
                obj.exp_info = rmfield(obj.exp_info,'flagged');
            catch
                disp('No flags found.')
            end
        end
        
        % Get the flagged tracks. 
        function flagged = get_flagged_tracks(obj)
            
            % Look for flagged tracks or cells. 
            if isfield(obj.exp_info.flagged,'tracks')
                flagged = obj.exp_info.flagged.tracks;
                warning('flagged tracks defined using old method!');
                
            elseif isfield(obj.exp_info.flagged,'cells')
                
                cells = obj.exp_info.flagged.cells;
                flagged = []
                %Loop over tracks. Look for flagged cells. 
                for i = 1:length(obj.tracks)
                    
                    for c= 1:size(cells,1)
                        
                        if obj.tracks{i}(1,1) == cells(c,1) & obj.tracks{i}(1,2) == cells(c,2)
                            flagged = [flagged, i];
                        end
                    end
                end
                
            else
                warning('no flagged objects')
                flagged=[];
                
            end
                
                

        end
        
        %quick plot of track length distribution. 
        function obj = distribution_track_length( obj, n )
            
            if nargin == 1
                n=20;
            end
            
            vals = cellfun(@(x) size(x,1), obj.tracks);
            figure(34);
            std_figure_settings
            hist(vals,n);
            
        end
        
        
        %Update the exp_info
        function obj = update_exp_info(obj, exp_info )
            
            %Check if combined experiment, regenerate paths to segmentation
            %files and image files. 
            if isfield(obj.exp_info,'combined_experiments')
                
                %Simple replacment. 
                names = setdiff( fieldnames(exp_info), {'seg_files','img_file','max_p_img'});
                for i = 1:length( names )
                    obj.exp_info = setfield(obj.exp_info,names{i},getfield(exp_info,names{i}));
                end
                
                %Seg files. 
                F = obj.exp_info.seg_files;
                for i = 1:length(F)
                    str_parts = strsplit(F{i},'/');
                    yes = strcmp(str_parts,'cell_lines');
                    start_idx = find(yes) + 1;
                    F{i} = fullfile(obj.exp_info.path_2_data,strjoin(str_parts(start_idx:end),'/'));
                end
                obj.exp_info.seg_files=F;               
                
                %Img files. 
                F = obj.exp_info.img_file;
                for i = 1:length(F)
                    str_parts = strsplit(F{i},'/');
                    yes = strcmp(str_parts,'cell_lines');
                    start_idx = find(yes) + 1;
                    F{i} = fullfile(obj.exp_info.path_2_data,strjoin(str_parts(start_idx:end),'/'));
                end
                obj.exp_info.img_file = F;
                
                %Maxp img files. 
                F = obj.exp_info.max_p_img;
                for i = 1:length(F)
                    str_parts = strsplit(F{i},'/');
                    yes = strcmp(str_parts,'cell_lines');
                    start_idx = find(yes) + 1;
                    F{i} = fullfile(obj.exp_info.path_2_data,strjoin(str_parts(start_idx:end),'/'));
                end
                obj.exp_info.max_p_img = F;
                
            else
                %Go through list of fields and update 
                names = fieldnames( exp_info );
                for i = 1:length( names )
                    obj.exp_info = setfield(obj.exp_info,names{i},getfield(exp_info,names{i}));
                end
            end
            
        end
        
        %Add image background data. 
        function obj = add_bg(obj,BG)
           
            %Check input. 
            if ~isstruct(BG)
                error('non-structure input');
            end
            
            %if isempty(obj.background)
            %    obj.background = BG;
            %else
            %    error('write code to merge background structures...')
            %end
            
        end
        
        %Remove data field from info. 
        function obj = rm_info_field(obj,field_str)
           
            obj.exp_info = rmfield(obj.exp_info,field_str);
        end
          
        % reset flaggs. 
        function obj = reset_flags(obj)
           
            obj = obj.rm_info_field('drawn_cells_frames');
            obj = obj.rm_info_field('drawn_cells');
            obj = obj.rm_info_field('flagged');
        end
        
        %Retrieve segmentation files
        function seg_files = get_frame_files(obj) 
            
            %First look in exp_info. If we combined experiments then there
            %should be a seg_files field in the exp_info. 
            if(isfield(obj.exp_info,'seg_files'))
                seg_files = obj.exp_info.seg_files;
            else
                f = dir([obj.exp_info.nuc_seg_dir,'/*.mat']);
                A = strvcat(f.name);
                B = repmat(obj.exp_info.nuc_seg_dir,[size(A,1),1]);
                seg_files = cellstr([B,A]);
            end
        end
        
        %Determine the default data channel. 
        function channel = findchannel(obj)
            %Check the structure.
            if isfield( obj.nuc_cyto_data, 'data')
                channel='data';
            elseif isfield( obj.nuc_cyto_data,'channel_01')
                channel='channel_01';
            else
                channel=[];
            end
        end
        
        %Get nuc_cyto data for a specific track. Can specify channel. 
        function data = get_track_data(obj, varargin)
            
            if length(varargin)==1
                idx = varargin{1};
                channel=obj.findchannel;
                track=[];
            else
                p = inputParser;
                p.addParameter('idx',[]);
                p.addParameter('channel',obj.findchannel,@isstring);
                p.addParameter('track',[]);
                p.parse(varargin{:});

                idx     = p.Results.idx;
                track   = p.Results.track;
                channel = p.Results.channel;
            end
            
            if ~isempty(idx)
                data_ids = obj.tracks{idx}(:,end);
            elseif ~isempty(track)
                data_ids = track(:,end);
            end
                      
            if isempty(channel)
                data = obj.nuc_cyto_data(data_ids);
            else
                data = obj.nuc_cyto_data.(channel)( data_ids );
            end

        end
        
        %Return all cell contours. 
        function cells = getCellContours(  obj, params )
            
            if nargin < 2
                params.seg_channel=1;
            end
            
            %Load frame obj. 
            seg_files = obj.get_frame_files;
            %Loop and load. 
            
            cells= cell(length(seg_files),1);
            disp_str = '';
            for i=1:length(seg_files)
                load( seg_files{i},'frame_obj');

                try
                    % Look for channel specific segmentation. 
                    ch_str = ['seg_channel_',pad(num2str(params.seg_channel),2,'left','0')];
                    if isfield(frame_obj,'contours')
                        
                        cells{i} = frame_obj.contours;
                    elseif isfield(frame_obj,ch_str)
                        cells{i} = frame_obj.(ch_str).contours;
                    else
                        error('Channel segmentation not defined');
                    end
                catch

                    if isfield(frame_obj,'PixelIdxList')             
                        %Trace boundaries. 
                        if length(size(frame_obj.BW))==3
                            BW = max(frame_obj.BW,[],3);
                            ctrs = cat(1,frame_obj.centroids{:});
                            ctrs = ctrs(:,1:2);
                        else
                            BW = frame_obj.BW;
                            ctrs = cat(1,frame_obj.centroids{:});
                        end
                        C = bwboundaries(BW);
                        if ~isempty(C)
                            C = cellfun(@(x) [smooth(x(:,2)),smooth(x(:,1))],C,'uniformoutput',0);
                            c_ctr = cellfun(@(x) [mean(x(:,1)),mean(x(:,2))],C,'uniformoutput',0);
                            %Match centroids. 
                            D = pdist2(cat(1,c_ctr{:}),ctrs);
                            [~,idx] = min(D);
                            frame_obj.contours=C( idx );
                            cells{i} = frame_obj.contours;
                            % Update frame_obj. 
                            save(seg_files{i},'frame_obj','-append');
                        end
                    else
                        cells{i} = [];
                    end
                end
                if exist(disp_str,'var'); clearString(disp_str); end
                disp_str=['loading frame ',num2str(i)];
                disp(disp_str);
            end          
        end
        
        %Get reader for image. 
        function reader = get_reader(obj)
           
            reader=bfGetReader(obj.exp_info.img_file);
            
        end
                
        % Generate an empty structure for keeping track of nuclear and cytoplasmic signals
        function data = gen_data_struct( obj, N )

            %Fields
            data(N) = struct('nuc_mean',[],'cyto_mean',[],'local_rho',[],'cell_id',[],'nuc_med',[],'cyto_med',[],'nuc_area',[]);

        end
        
        %Return the associated spot_tracks for a cell track. 
        function ids = get_spot_track_ids( obj, track_idx)
           
            %All cell track idxs 
            all_ids = cellfun(@(x) x(1,3), obj.spot_tracks);
            
            ids = find( all_ids == track_idx);
            
        end
        
        %Return intensity traces of spot tracks
        function int_traces = getSpotTrackIntTraces( obj, indices, int_thresh )

            %Loop over indices
            int_traces=cell(length(indices),1);
            indices= indices(:)';
            c=0;
            for i = indices

                %Frames. 
                frames = obj.spot_tracks{i}(:,1);

                %fit result ids
                result_ids = obj.spot_tracks{i}(:,2);

                %Intensities. 
                int = cat(1,obj.results( result_ids ).sum_int);
                c=c+1;
                int_traces{c} = [frames, int];
            end
            
            %If intensity threshold supplied. 
            c=0;
            filt_traces={};
            if nargin==3
                for i = 1:length(int_traces)
                    int = int_traces{i}(:,2);
                    sel = int >= int_thresh;
                    %Skip track if empty. 
                    if(sum(sel)==0)
                        continue
                    end
                    c=c+1;
                    filt_traces{c} = int_traces{i}(sel,:);
                end
                int_traces = filt_traces;
            end
                
        end
       
    end
    
end

