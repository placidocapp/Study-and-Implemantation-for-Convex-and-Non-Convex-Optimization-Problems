clear all;
close all;
clc

format long;
%% Inicialização

%Backtracking
ro = 0.9;
c = 10^-4;

%Todos
maxIter = 100;

%% Algoritmo backtraking

% syms x;
% f = f(x);
% df = diff(f);

%Inicializa
x = zeros(maxIter,1);

x(1) = 10;
for k = 1:maxIter
    d = -df(x(k));
    alpha = 1;
    while f(x(k) + alpha*d) > f(x(k)) + c*alpha*df(x(k))*d 
        alpha = ro*alpha;
    end
    x(k+1) = x(k) + alpha*d;
end
x
sol = f(x(end))


%% Zomm Algorithm with Strong Wolfe conditions

%Inicialization of x
x = zeros(maxIter,1);
x(1) = 0.5;

%Main loop
for k = 1:maxIter
    d = -df(x(k))  %Choose the direction as an steepest gradient descent

   alpha = strongwolfe(f,d,0.5,10);
    
    x(k+1) = x(k) + alpha*d;
    x(k+1)
end
    
%Display the x history and the solution value
x
sol = f(x(end))


