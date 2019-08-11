function [ alpha ] = backtrack( f, g, x, d, c1, ro )
%This function use the backtracking to find the step lenght

alpha = 1;      %Step size, firs gess
maxIter = 100;

k = 1;
%Reduce the step until it satisfies the Armijo condition
while f(x + alpha*d) > f(x) + c1*alpha*g(x)'*d 
    alpha = ro*alpha;
    k = k+1;
    if k == 100
        return
    end
end

end

