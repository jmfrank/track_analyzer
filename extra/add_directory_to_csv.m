function add_directory_to_csv()
%% Adding a directory of files to CSV

loc = 'Work'


csv_files{1} = '/Users/mf2741/Dropbox/TEAD_paper/stowers.csv';
csv_files{2} = '/home/matt/Dropbox/TEAD_paper/stanford_experiments.csv';
csv_files{3}  = '/Users/matt/Dropbox/TEAD_paper/stanford_experiments.csv';
csv_files{4} =  '/home/matt/Dropbox/TEAD_paper/stowers.csv';
csv_file=csv_files{2}


%Read csv file
[fid message ] = fopen(strtrim(csv_file));
data = textscan(fid,'%s %s',800,'delimiter',',');
fclose(fid);
start_row = size(data{1},1)+1;

%Get experiment set base dir
if( loc == 'Work')
    tmp = cellfun(@(x) x(3), data(1));
elseif( loc == 'Home')
    tmp = cellfun(@(x) x(4), data(1));
else
    disp('Location must be "Home" or "Work"');
end

selpath = uigetdir('/media/matt/internal_8TB/data/microscopy/stanford/');

% Loop through all subdirectories, if any. Dir is recursive using **. Just
% limit to ns2 file for now. 
file_types={'.tif','.nd2','.czi'}; %add more as necessary

all_files = get_files_of_type(selpath,file_types);
all_files = all_files(~[all_files.isdir]);
all_files
%Counter
c=1;

% pre-fix for the directory:
prefix=tmp{1};
for i = 1:length(all_files)
    %Add a row for this file
    data{1}{start_row+c} = erase(all_files(i).folder, prefix);
    data{2}{start_row+c} = all_files(i).name;
    c = c+ 1;
end   


%% First backup old csv
[a,b,c] = fileparts(csv_file);
copyfile(csv_file,fullfile(a,[b,'_old.csv']))
%% Update csv by converting data to cell matrix
new_data = cell( size(data{1},1),size(data,2));
for i = 1:length(data)
    L = size(data{i},1);
    new_data(1:L,i) = data{i};
end
%Print csv
cell2csv(csv_file,new_data,',');

end

function all_files = get_files_of_type(selpath,file_types)

    all_files=[];
    for i=1:length(file_types)
        these_files=dir2(fullfile(selpath, ['**/*',file_types{i}]));
        all_files=[all_files; these_files];
    end

end


