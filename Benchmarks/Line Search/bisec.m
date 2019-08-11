function [ t ] = bisec( f, g, x, p, c1, c2 )
%This functions chooses the step lenght that satisfies the weak wolfe
%conditions
% Receives the function f and it's gradient g. And the point x and
% direction p
% This algorithm was found in a lecture from Professor James V. Burke from
% University of Washington, Seattle. Here is the link for the algorithm 
%(page 8): https://sites.math.washington.edu/~burke/crs/408/notes/nlp/line.pdf
%and the autors page https://sites.math.washington.edu/~burke/ . This links
%are working on 07/10/2018
maxSubIter = 30; %Max iterations of the bisection algorithm

beta = -1;
alpha = 0;
t = 1;
i = 0;
gk = g(x);
while(1)
    if f(x+t*p) > f(x) + c1*t*gk'*p
        dalpha = gk'*p;
        t = (-dalpha*beta^2)/( 2*(f(x+t*p)-f(x)-dalpha*beta) );
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

