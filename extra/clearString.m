%Remove string from displace so we can write a new string. 
function clearString(disp_str)

    fprintf(repmat('\b',length(disp_str)+1));

end