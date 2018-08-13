function [ T ] = changebase( T, i, j )
%Insert a new base in position (i,j)

T(i,:) = T(i,:)/T(i,j);
for k = 1:size(T,1)
    if k ~= i
        T(k,:) = T(k,:) - T(i,:)*(T(k,j));
    end
end

end

