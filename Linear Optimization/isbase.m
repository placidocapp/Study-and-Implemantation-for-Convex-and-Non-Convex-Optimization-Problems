function [ y ] = isbase( v )
%Verify if the column is a base, return 0 or 1

if (sum(v == 0) == length(v)-1) && (sum(v == 1) == 1)
    y = 1;
else
    y = 0;
end

end

