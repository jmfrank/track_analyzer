function bw = border_frame( S )

bw = zeros(S);

if length(S) == 2
       
    bw([1,end],:) = 1;
    bw(:,[1,end]) = 1;     
    
elseif length(S) == 3
    
    bw([1,end],:,:) = 1;
    bw(:,[1,end],:) = 1;
    %bw(:,:,[1,end]) = 1;
    
else
    
    error('Not supported input.')
    
end

end
