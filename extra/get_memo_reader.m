% Retrieve data obj, or generate one if it doesn't exist. 
% info structure must contain csv_file and exp_id (row in csv_file)

function reader = get_memo_reader(img_file) 

    %Creater reader.
    reader = loci.formats.Memoizer(bfGetReader(),0);
    reader.setId(img_file);


end
