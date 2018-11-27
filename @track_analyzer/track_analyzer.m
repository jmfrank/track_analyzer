classdef track_analyzer

    
    properties (SetAccess = private)
        
        %Store the experimental info (i.e. file directories, experimental parameters, etc. 
        exp_info = {};
        
        % The trajectories stored in a cell array, one T x n_dim per particle
        tracks = {};
        % Trajectories corresponding to the spots localized in cell tracks
        spot_tracks = {};
        
        %Cells (originally, 2D contour traces). Need to change for 3D
        %segmentation. 
        cells
        
        %Background value. Each experiment should have a background value
        %in order to calculate intensity properly. 
        img_bg 
        
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
    end
    
    %Hidden variables for checking things
    properties (SetAccess = private, Hidden = true)
        nuc_cyto_calc = false;
        nuc_area_calc = false;
        seg           = false;
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
        
        %Remove data field from info. 
        function obj = rm_info_field(obj,field_str)
           
            obj.exp_info = rmfield(obj.exp_info,field_str);
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
                p.addParameter('track',[],@iscell);
                p.parse(varargin{:});

                idx = p.Results.idx;
                track = p.Results.track;
                channel = p.Results.channel;
            end
            
            
            if ~isempty(idx)
                data_ids = obj.tracks{idx}(:,end);
            elseif ~isempty(track)
                data_ids = track{1}(:,end);
            end
            
            if isempty(channel)
                data = obj.nuc_cyto_data(data_ids);
            else
                data = obj.nuc_cyto_data.(channel)( data_ids );
            end


        end
        
        %Return all cell contours. 
        function cells = getCellContours(  obj )
            %Load frame obj. 
            seg_files = obj.get_frame_files;
            %Loop and load. 
            
            cells= cell(length(seg_files),1);
            disp_str = '';
            for i=1:length(seg_files)
                tmp = load( seg_files{i},'frame_obj');

                try

                    cells{i} = tmp.frame_obj.contours;
                catch

                    if isfield(tmp.frame_obj,'PixelIdxList')             
                        cells{i} = tmp.frame_obj.PixelIdxList;
                    else
                        cells{i} = [];
                    end
                end
                if exist(disp_str,'var'); clearString(disp_str); end;
                disp_str=['loading frame ',num2str(i)];
                disp(disp_str);
            end          
        end
        
        % Generate an empty structure for keeping track of nuclear and cytoplasmic signals
        function data = gen_data_struct( N )

            %Fields
            data(N) = struct('nuc_mean',[],'cyto_mean',[],'local_rho',[],'cell_id',[],'nuc_med',[],'cyto_med',[],'nuc_area',[]);

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

