Algorithm from Numerical Optimization. Jorge Nocedal and Stephen J. Wright.
Second edition. Zoom algorithm can be found at pages 60 and 61.
MATLAB code by Plácido Campos
Last modified: may 28, 20018

clear all;
close all;
clc

format long;
%% Parameters 

% Wolfe conditions
c1 = 10^-4;
c2 = 0.5;
alphaMax = 10;
eps = 10^-20;

%Backtracking
ro = 0.8;

%Todos
maxIter = 100;

%% Backtraking Algorithm

%Use this to find de derivative of a new function
% syms x;
% f = f(x);
% df = diff(f);

%Inicialization of x
x = zeros(maxIter,1);
x(1) = 0.5;

%Backtraking loop
for k = 1:maxIter
    d = -df(x(k));  %Choose the direction as an steepest gradient descent
    alpha = 1;      %Step size, firs gess
    
    %Reduce the step until it satisfies the Armijo condition
    while f(x(k) + alpha*d) > f(x(k)) + c1*alpha*df(x(k))*d 
        alpha = ro*alpha;
    end
    
    %Find the next point
    x(k+1) = x(k) + alpha*d;
end

%Display the x history and the solution value
x
sol = f(x(end))
disp('These are te results of the backtraking, press anything to zoom...')
pause

%% Zomm Algorithm with Strong Wolfe conditions

%Inicialization of x
x = zeros(maxIter,1);
x(1) = 0.5;

%Main loop
for k = 1:maxIter
    d = -df(x(k))  %Choose the direction as an steepest gradient descent

    %Inicializations for alpha search
    i = 1;
    alpha = 1;      %Initial gess for alpha
    alpha_ant = 0;
    while (1)
        fi = f(x(k) + d*alpha);     %Contains the fi(alpha) value
        fio = f(x(k));              %fi(0)
        dfio = df(x(k));            %dfi(0)
        if (fi + eps*sign(fi) > fio + c1*alpha*dfio)||(i>1 && fi >= f(x(k) + d*alpha_ant) )
            disp('Condição 1')
            alpha = zoom(alpha_ant, alpha, x(k), d, c1, c2, eps);
            break;
        end
        dfi = df(x(k) + d*alpha)   %Contains the dfi(alpha) value
        if abs(dfi)+ eps >= -c2*dfio
            disp('Condição 2')
            break;
        end
        if dfi >= 0
            disp('Condição 3')
            alpha = zoom(alpha, alpha_ant, x(k), d, c1, c2, eps);
            break;
        end
        
        alpha_ant = alpha;
        alpha = (alphaMax+alpha)/2;
        i = i + 1;
    end
    
    x(k+1) = x(k) + alpha*d;
    x(k+1)
end
    
%Display the x history and the solution value
x
sol = f(x(end))    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



