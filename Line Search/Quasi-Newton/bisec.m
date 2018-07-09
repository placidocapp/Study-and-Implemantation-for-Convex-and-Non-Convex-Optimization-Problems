function [ t ] = bisec( f, g, x, p, c1, c2 )
%This functions chooses the step lenght that satisfies the weak wolfe
%conditions
% Receives the function f and it's gradient g. And the point x and
% direction p

maxSubIter = 30; %Max iterations of the bisection algorithm

beta = -1;
alpha = 0;
t = 1;
i = 0;
gk = g(x);
while(1)
    if f(x+t*p) > f(x) + c1*t*gk'*p
        t = -g(x)*beta^2/( 2*(f(x+t*p)-f(x)-g(x)*beta) );
    elseif g(x+t*p) < c2*gk'*p
        alpha = t;
        if beta == -1
            t = 2*alpha;
        else
            t = 0.5*(alpha+beta);
        end
    else
        break;
    end
    %Stop creteria to avoid being stuck in case the algorithm choosed an
    %increase direction
    i = i+1;
    if i >= maxSubIter
        t = 0.001;
        break;
    end
end

end

