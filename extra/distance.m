function C = distance(A,B)

% calculate pairwise linear distance between each points. 

C = inf(size(A));

for i = 1:size(A,1)

    for j = 1:size(B,1)

        a=A(i,:);
        b=B(j,:);
        if length(a) == 2

            C(i,j) = (a(1)-b(1))^2 + (a(2)-b(2))^2 ; %+ (a(3)-b(3))^2;

        end

    end
end
