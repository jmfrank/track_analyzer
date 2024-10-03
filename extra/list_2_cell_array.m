% covert a stats structure list to a cell array. 
function out_cell = list_2_cell_array(L,field_name)
    out_cell=cell(length(L),1);
    for i = 1:length(L)       
        out_cell{i} = L(i).(field_name);
    end

end