function [ t ] = zoom( g, f, x, p, l, h, c1, c2, maxSubIter )
%Receives a low and a high point and zoom between this values to find the
%point that satisfies the strong wolfe conditions

i = 1;

while(1)
    %% Initial gess using bisection
    
    t = 0.5*(l+h);
    
    %% Calculate the function values
    fak = f(x+t*p);
    fk = f(x);
    gak = g(x+t*p)'*p;
    gk = g(x)'*p;

    if fak > fk+c1*t*gk || fak >= f(x+l*p)
        h = t;
    else
        if abs(gak) <= -c2*gk
            return;
        end
        
        if gak*(h-l) >= 0
            h = l;
        end
        l = t;
    end
    
   
    
    %% Safeguard Stop Criteria
    
    i = i + 1;
    if i > maxSubIter
        break
    end
end

