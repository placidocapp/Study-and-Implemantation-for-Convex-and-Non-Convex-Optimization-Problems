function [ alpha ] = zoom( alphalow, alphahi, x, d, c1, c2, eps )
%This function will zoom over the area were the biggest step that satisfies
%the wolf conditions must be  
%   Receive the alphas high and low (the limits of best alpha) and the
%   actual x
disp('entrou fun��o')
fio = f(x)             %fi(0)
dfio = df(x)           %dfi(0)
while (1)
    %Interpolate using a quadratic function to fint a point between alphahi 
    %alphalow. (We choose the minimum point of the quadratic)
    %alpha = alphalow - df(x+alphalow*d)*(alphahi-alphalow)^2/...
    %    ( 2*( f(x+d*(alphahi-alphalow)) -f(x+alphalow*d)...
    %   -df(x+alphalow*d)*(alphahi-alphalow) ) )
    
    %TESTE
    alphalow
    alphahi
    
    %Interpolate using bisection approach
    alpha = (alphahi + alphalow)/2
    
    fi = f(x + d*alpha)    %fi(alpha)
    
    
    %TESTE
    fi_low = f(x + d*alphalow)
    fi_hi = f(x + d*alphahi)
    
    %This part was added by me, sometimes the algorithm tends to the value
    %that satisfies the conditions but never reaches that point, so if its
    %near enough (difference bigger than eps) than stop
%     if abs(alphahi-alphalow) < eps
%         alpha = (alphahi+alphalow)/2;
%         break;
%     end

    disp('condi��es')
    c1 = fi
    c2 = fio + c1*alpha*dfio
    c3 = fi - eps*sign(fi)
    c4 = f(x + d*alphalow)
    if (fi + eps*sign(fi) > fio + c1*alpha*dfio)||(fi + eps*sign(fi) >= f(x + d*alphalow))  
        disp('Condi��o Zoom 1')
        alphahi = alpha;
    else
        dfi = df(x + d*alpha)
        -c2*dfio
        if abs(dfi) + eps >= -c2*dfio
            disp('Condi��o Zoom 2')
            break;
        elseif dfi*(alphahi - alphalow) >= 0
            disp('Condi��o Zoom 3')
            alphahi = alphalow;
        end
        disp('passou direto')
        alphalow = alpha;
    end
end
disp('Saiu fun��o')

end

