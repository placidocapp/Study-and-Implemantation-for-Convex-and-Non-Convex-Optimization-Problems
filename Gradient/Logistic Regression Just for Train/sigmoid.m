function [ y ] = sigmoid( z )
%Calculate the sigmoid function

l = size(z,1);      %Number of lines
c = size(z,2);      %Number of Columns
y = zeros(l,c);     %Size of y is the same as z

%Calculate the sigmoid function for each term
for i = 1:l
    for j = 1:c
        y(i,j) = 1/(1+exp(-z(i,j)));
    end
end

end

