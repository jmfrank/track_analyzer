%Retrieve info on experiment specified in a csv_file. 
%This was written for my own specific organization style. Adapt as needed. 
% exp_info output needs to be passed to track_analyzer to initialize data tracking object. 

function exp_info = get_exp_info(info) 

%Need to figure out which computer we are using. I.e. which path to
%dropbox. 
try
    full_file = fullfile('/Users/franklin/Dropbox/',info.csv_file);
    %Read csv file
    [fid message ] = fopen(strtrim(full_file));
    data = textscan(fid,'%s %s %s %s %s %s %q %s %s',800,'delimiter',',');
    fclose(fid);
    loc = 'Home';
catch
    full_file = fullfile('/home/jan/Dropbox',info.csv_file);
     %Read csv file
    [fid message ] = fopen(strtrim(full_file));
    data = textscan(fid,'%s %s %s %s %s %s %q %s %s',800,'delimiter',',');
    fclose(fid);
    loc = 'Work';
end

%If specified, use the input location
if isfield(info,'loc')
    loc = info.loc;
end


%Get experiment set base dir
if strcmp(loc , 'Work')
    tmp = cellfun(@(x) x(3), data(1));
elseif strcmp( loc , 'Home')
    tmp = cellfun(@(x) x(4), data(1));
elseif strcmp( loc, 'local')
    tmp = cellfun(@(x) x(5), data(1));
elseif strcmp( loc,'dropbox')
    a = cellfun(@(x) x(6), data(1));
    tmp{1} = [path_2_dropbox,a{1}];
else
    
    disp('Location must be "Home" or "Work"');
end
%now fill in extra columns of data if there's a mis-match at the end. 
L = cellfun(@(x) size(x,1),data);
max_L = max(L);
bad_columns = find(L < max_L);
for i = 1:length(bad_columns)
    data{bad_columns(i)}{max_L} = {''};
end

base_dir = tmp{1};
%Get specified experiment info
all_info = cellfun(@(x) x(info.exp_id), data);

%Specify variables <<<<*** ADD MORE LATER??? *** >>>>
exp_info.sub_dir     = strcat(base_dir,all_info{1});


%Use the 'type' field to figure out how to parse info. 
if ~isfield( info, 'type' )
    data_type='standard';
else
    data_type= info.type;
end


switch data_type
    
    
    case 'standard'
        

        %check if img_file is empty
        if( isempty( all_info{2} ) )

            %Now ask user which file to select in the sub_dir
            exp_info.sub_dir
            [sel_file,sel_path,indx] = uigetfile(fullfile(exp_info.sub_dir,'*.*'));

            %Create the img_file name
            exp_info.img_file = [sel_path,'/', sel_file];

            %Now we need to update the csv file?
            data{ 2 }{exp_id} = sel_file;

            %Convert data to cell matrix
            new_data = cell( size(data{1},1),size(data,2));
            for i = 1:length(data)
                L = size(data{i},1);
                new_data(1:L,i) = data{i};
            end
            %Print csv
            cell2csv(csv_file,new_data,',');

        else
            exp_info.img_file    = fullfile(base_dir,all_info{1},all_info{2});

        end
    
    
        [a,b,c] = fileparts( exp_info.img_file );

        exp_info.track_file  = fullfile(a,[b,'_track_data.mat']);

        %Make a directory based on the lsm_file name.  
        [a,b,c] = fileparts( exp_info.img_file );
        seg_dir = fullfile(a,[b,'/']);
        if( ~exist( seg_dir, 'dir') )
            mkdir( seg_dir );
        end
         exp_info.nuc_seg_dir = seg_dir;

         %Max projection file. 
         exp_info.max_p_img = fullfile(a,[b,'_maxp.tif']);
    
    %user specified interpretation. 
    case 'pre_img_frap'
        
        exp_info.pre_img_file    = fullfile(base_dir,all_info{1},all_info{2});
        [a,b,c] = fileparts( exp_info.pre_img_file );
        seg_dir = fullfile(a,[b,'/']);
        if( ~exist( seg_dir, 'dir') )
            mkdir( seg_dir );
        end
        exp_info.pre_img_nuc_seg_dir = seg_dir;
        exp_info.max_p_pre_img = fullfile(a,[b,'_maxp.tif']);

        
        exp_info.img_file    = fullfile(base_dir,all_info{1},all_info{3});
        [a,b,c] = fileparts( exp_info.img_file );
        seg_dir = fullfile(a,[b,'/']);
        if( ~exist( seg_dir, 'dir') )
            mkdir( seg_dir );
        end
        exp_info.nuc_seg_dir = seg_dir;
        %Max projection file. 
        exp_info.max_p_img = fullfile(a,[b,'_maxp.tif']);
        
        %track
        exp_info.track_file  = fullfile(a,[b,'_track_data.mat']);

end

    
%Add path to data. 
exp_info.path_2_data = base_dir;

end
