function [ alpha ] = zoom( alphalow, alphahi, x, d, c1, c2, eps )
%This function will zoom over the area were the biggest step that satisfies
%the wolf conditions must be  
%   Receive the alphas high and low (the limits of best alpha) and the
%   actual x

disp('entrou função')
while (1)
    %Interpolate using a quadratic function to fint a point between alphahi 
    %alphalow. (We choose the minimum point of the quadratic)
    %alpha = alphalow - df(x+alphalow*d)*(alphahi-alphalow)^2/...
    %    ( 2*( f(x+d*(alphahi-alphalow)) -f(x+alphalow*d)...
    %   -df(x+alphalow*d)*(alphahi-alphalow) ) )
    
    %TESTE
    alphalow;
    alphahi;
    
    %Interpolate using bisection approach
    alpha = (alphahi + alphalow)/2;
    
    fi = f(x + d*alpha);    %fi(alpha)
    fio = f(x);             %fi(0)
    dfio = df(x);           %dfi(0)
    
    %TESTE
    fi_low = f(x + d*alphalow);
    fi_hi = f(x + d*alphahi);
    
    %This part was added by me, sometimes the algorithm tends to the value
    %that satisfies the conditions but never reaches that point, so if its
    %near enough (difference bigger than eps) than stop
    if abs(alphahi-alphalow) < eps
        alpha = (alphahi+alphalow)/2;
        break;
    end
    if (fi + eps*sign(fi) > fio + c1*alpha*dfio)||(fi + eps*sign(fi) >= f(x + d*alphalow))  
        disp('Condição Zoom 1')
        alphahi = alpha;
    else
        dfi = df(x + d*alpha);
        c2*dfio*d;
        if abs(dfi)+ eps >= -c2*dfio
            disp('Condição Zoom 2')
            break;
        elseif dfi*(alphahi - alphalow) >= 0
            disp('Condição Zoom 3')
            alphahi = alphalow;
        end
        disp('passou direto')
        alphalow = alpha;
    end
    
end
disp('Saiu função')

end

