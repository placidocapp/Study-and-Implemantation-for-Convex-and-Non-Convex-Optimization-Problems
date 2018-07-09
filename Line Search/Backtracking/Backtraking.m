% Algorithm from Numerical Optimization. Jorge Nocedal and Stephen J. Wright.
% Second edition. Zoom algorithm can be found at pages 60 and 61.
% MATLAB code by Plácido Campos
% Last modified: may 30, 20018

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
ro = 0.5;

%Todos
maxIter = 100;

%Size
n = size(df(zeros(10,10)),1);

%% Backtraking Algorithm

%Use this to find de derivative of a new function
% syms x;
% f = f(x);
% df = diff(f);

%Inicialization of x
x = zeros(maxIter,n);
x(1,:) = randn(1,n);

%Backtraking loop
for k = 1:maxIter
    d = -df(x(k,:));  %Choose the direction as an steepest gradient descent
    alpha = 1;      %Step size, firs gess
    
    %Reduce the step until it satisfies the Armijo condition
    while f(x(k,:)' + alpha*d) > f(x(k,:)') + c1*alpha*df(x(k,:))'*d 
        alpha = ro*alpha;
    end
    
    %Find the next point
    x(k+1,:) = (x(k,:)' + alpha*d)';
end

%Display the x history and the solution value
x
sol = f(x(end,:)')