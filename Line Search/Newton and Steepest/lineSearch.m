function [ t ] = lineSearch( f, g, x, p, c1, c2 )
%Call zoom algorithm

%% Initializations

alphaMax = 10;      %The superior limit
i = 2;              %Counter of iterations
maxSubIter = 30;    %Maximum iterations for this algorithm
alpha = zeros(maxSubIter+1,1);
alpha(1) = 0;       %The inferior limit
alpha(2) = 1;       %The first iteration

%% Algorithm

while(1)
    %% Calculate all values of the functions
    
    fak = f(x+alpha(i)*p);
    fk = f(x);
    gak = g(x+alpha(i)*p)'*p;
    gk = g(x)'*p;
 
    %% Chose to call or not the zoom
    if fak > fk+c1*alpha(i)*gk || (fak >= f(x+alpha(i-1)*p) && i > 2)
        t = zoom(g,f,x,p,alpha(i-1),alpha(i),c1,c2,maxSubIter);
        return;
    end
    
    if abs(fak) <= -c2*gk
        t = alpha(i);
        return;
    end
    
    if gak >= 0
        t = zoom(g,f,x,p,alpha(i),alphaMax,c1,c2,maxSubIter);
        return;
    end
    
    %% Update
    
    alpha(i+1) = 0.5*( alpha(i) + alphaMax );
    i = i + 1;

    %% Safeguard Stop Criteria
    if i > maxSubIter
        t = alpha(i);
        break
    end
end

