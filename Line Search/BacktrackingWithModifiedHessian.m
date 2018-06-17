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

%Size (Selects size based on the gradient size)
n = size(df(zeros(10,10)),1);

%% Backtraking Algorithm (Modified Hessian with 2-norm)

%Use this to find de derivative of a new function
% syms x;
% f = f(x);
% df = diff(f);

%Inicialization of x
x = zeros(maxIter,n);
x(1,:) = randn(1,n);
% x(1,:) = 0.333510833065806;

%Backtraking loop
for k = 1:maxIter
    
    %Calculates the hessian, if some eigenvalue of H is negative than
    %change the hessian for itself plus eye*|most negative eigenvalue + lambda|
    H = d2f(x(k,:));
    if sum(eig(H) < 0) > 0
        lambda = min(eig(H));
        H = H + eye*(-lambda+0.1);
    end
    
    d = -inv(H)*df(x(k,:));  %Choose the direction as an steepest gradient descent
    alpha = 1;      %Step size, firs gess
    
    %Reduce the step until it satisfies the Armijo condition
    while f(x(k,:)' + alpha*d) > f(x(k,:)') + c1*alpha*df(x(k,:))'*d 
        alpha = ro*alpha;
    end
    
    %Find the next point
    x(k+1,:) = (x(k,:)' + alpha*d)';
end
out = [x, x, x];
sol = [f(x(end,:)') f(x(end,:))' f(x(end,:))'];

%% Backtraking Algorithm (Modified Hessian with cholesky factorization)

%Use this to find de derivative of a new function
% syms x;
% f = f(x);
% df = diff(f);

%Inicialization of x
x = zeros(maxIter,n);
x(1,:) = randn(1,n);
% x(1,:) = 0.333510833065806;

%Backtraking loop
for k = 1:maxIter
    
    %Calculates the hessian, if some eigenvalue of H is negative than
    %change the hessian for itself plus eye*|most negative eigenvalue + lambda|
    H = d2f(x(k,:));
    if sum(eig(H) < 0) > 0
        [L,D] = mcfac(H);
        H = L*D*L';
    end
    
    d = -inv(H)*df(x(k,:));  %Choose the direction as an steepest gradient descent
    alpha = 1;      %Step size, firs gess
    
    %Reduce the step until it satisfies the Armijo condition
    while f(x(k,:)' + alpha*d) > f(x(k,:)') + c1*alpha*df(x(k,:))'*d 
        alpha = ro*alpha;
    end
    
    %Find the next point
    x(k+1,:) = (x(k,:)' + alpha*d)';
end
out(:,2) = x
sol(:,2) = f(x(end,:)')
%% Backtraking Algorithm (Steepest Descent)

%Use this to find de derivative of a new function
% syms x;
% f = f(x);
% df = diff(f);

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
    
disp('Compare the steepest gradient descent vs newton method:')
out(:,3) = x
sol(:,3) = f(x(end,:)')
    
%Plot the error
error = zeros(maxIter,2);
for j = 1:3
    for i = 1:maxIter
        error(i,j) = abs(f(out(i,j)) - mean(sol));
    end
end

plot(error(:,1), 'b'), hold on, plot(error(:,2), 'g'), plot(error(:,3), 'r')
legend('Newton Method With 2-norm','Newton Method With Cholesky','Steepest Gradient Descent')
title('Error of Steepest Gradient Descent vs Newton Method')
% axis([0 10 0 5])
    
    
    
    
    
    
    
    
    
    
    
    
    
    



