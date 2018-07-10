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
eps = 10^-8;

%Backtracking
ro = 0.9;

%Todos
maxIter = 200;

%Size
n = 2;

kfinal = -1;

%% Function

%Function values
A = randn(n);
A = A'*A;
A = A'+ A;
b = randn(n,1);

%Calculate the function value and its derivatives
f = @(x) 0.5*x'*A*x + b'*x;
g = @(x) A*x + b;
B = @(x) A;
sol = - (pinv(A)*b)';

% f = @(x) sin(x(1))+sin(x(2));
% g = @(x) [cos(x(1));
%           cos(x(2))];
% B = @(x) [  -sin(x(1)) 0;
%             0          -sin(x(2))];
% sol = [-1.570796327268957  -1.570796324699814];

% f = @(x) -200*exp(-0.2*sqrt(x(1)^2+x(2)^2));
% g = @(x) [(4*x(1).*exp(-(x(1).^2 + x(2).^2)^(1/2)/50))./(x(1).^2 + x(2).^2)^(1/2);
%           (4*x(2).*exp(-(x(1).^2 + x(2).^2)^(1/2)/50))./(x(1).^2 + x(2).^2)^(1/2)];
% B = @(x) [ (4*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(1/2) - (2*x(1)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2),                                                    - (2*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2);
%                                                     - (2*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2), (4*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(1/2) - (2*x(2)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(2)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2)] ;
% sol = [0 0];

%Rosenbrock
% a = 1;
% b = 100;
% f = @(x) (a-x(1))^2 + b*(x(2)-x(1)^2)^2; 
% g = @(x) [-2*(a - x(1))-4*b*x(1)*(x(2)-x(1)^2);
%           2*b*(x(2)-x(1)^2)];
% B = @(x) [2-4*b*x(2)+12*b*x(1)^2 -4*b*x(1);
%           -4*b*x(1)                2*b];       
% sol = [a a^2];

%% Backtraking Algorithm

%Inicialization of x
x = zeros(maxIter,n);
x(1,:) = randn(1,n);

%Backtraking loop
for k = 1:maxIter
    d = -g(x(k,:)');  %Choose the direction as an steepest gradient descent
    alpha = 1;      %Step size, firs gess
    
    %Reduce the step until it satisfies the Armijo condition
    while f(x(k,:)' + alpha*d) > f(x(k,:)') + c1*alpha*g(x(k,:)')'*d 
        alpha = ro*alpha;
    end
    
    %Find the next point
    x(k+1,:) = (x(k,:)' + alpha*d)';
    
    % Stop criteria
    if norm(g(x(k,:)')) < eps
        disp('Stopped at')
        kfinal = k
        break
    end
    
end

if kfinal == -1
    kfinal = maxIter;
end

x(kfinal,:)
norm(g(x(kfinal,:)'))
f(x(kfinal,:)')

%% Plot some graphs

x_error = zeros(kfinal,1);
f_error = zeros(kfinal,1);
g_evol = zeros(kfinal,1);
for i = 1:kfinal
    x_error(i) = norm(sol) - norm(x(i,:));
    f_error(i) = f(sol') - f(x(i,:)');
    g_evol(i) = norm(g(x(i,:)'));
end

k = 1:kfinal;

subplot(3,1,1), plot(k,x_error), title('Error in position x')
subplot(3,1,2), plot(k,f_error), title('Error in function value')
subplot(3,1,3), plot(k,g_evol), title('Gradient evolution')
    

