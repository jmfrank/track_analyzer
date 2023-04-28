
%% Adding a directory of files to CSV
clear

loc = 'Work'

%csv_file = '/home/matt/Dropbox/TEAD_paper/data/tead_experiments.csv'
csv_file = '/Users/mf2741/Dropbox/TEAD_paper/stowers.csv';

%Read csv file
[fid message ] = fopen(strtrim(csv_file));
data = textscan(fid,'%s %s %s %s %s %s %s',800,'delimiter',',');
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

selpath = uigetdir('/media/matt/DATA_1_guan/cell_lines');

% Loop through all subdirectories, if any. Dir is recursive using **. Just
% limit to ns2 file for now. 
these_files = dir2(fullfile(selpath, '**/*.czi'));
these_files = these_files(~[these_files.isdir]);

%Counter
c=1;

% pre-fix for the directory:
prefix=tmp{1};
for i = 1:length(these_files)

    %Add a row for this file
    data{1}{start_row+c} = cut_string_at(these_files(i).folder);
    data{2}{start_row+c} = these_files(i).name;
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

