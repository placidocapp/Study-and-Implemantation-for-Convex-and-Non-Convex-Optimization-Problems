function [ alpha ] = backtrack( f, g, x, d, c1, ro )
%This function use the backtracking to find the step lenght

 alpha = 1;      %Step size, firs gess
%Reduce the step until it satisfies the Armijo condition
while f(x + alpha*d) > f(x) + c1*alpha*g(x)'*d 
    alpha = ro*alpha;
end

end

