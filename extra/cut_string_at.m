

%% Cut a path after a specific string

function cut_str = cut_string_at(in_str, match)

str_parts = strsplit(in_str,'/');
yes = strcmp(str_parts,'cell_lines');
start_idx=find(yes)+1;
cut_str = strjoin(str_parts(start_idx:end),'/');

end