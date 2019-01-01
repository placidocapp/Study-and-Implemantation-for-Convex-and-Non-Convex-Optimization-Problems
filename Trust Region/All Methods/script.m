% MATLAB code by Plácido Campos based on Jorge Nocedal and Stephen J. Wright.
% Numerical Optimization, Second Edition. 

clear all;
close all;
clc

format long;
%% Parameters 

%Stop criteria for norm of gradient
eps = 10^-8;

%Trust region constants
eta = 0.25;
r = 1;

n = 2;

%Limit of iterations
maxIter = 100;

%% Function

%Function values
A = randn(n);
A = A'*A;
A = A'+ A;
b = randn(n,1);

%Calculate the function value and its derivatives
% f = @(x) 0.5*x'*A*x + b'*x;
% g = @(x) A*x + b;
% B = @(x) A;
% sol = - (pinv(A)*b)';

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
a = 1;
b = 100;
f = @(x) (a-x(1))^2 + b*(x(2)-x(1)^2)^2; 
g = @(x) [-2*(a - x(1))-4*b*x(1)*(x(2)-x(1)^2);
          2*b*(x(2)-x(1)^2)];
B = @(x) [2-4*b*x(2)+12*b*x(1)^2 -4*b*x(1);
          -4*b*x(1)                2*b];       
sol = [a a^2];

%% Chama algoritmo

opt = options('cauchy','bfgs');
%algorithm     = 'cauchy', 'dogleg', 'subproblem'
%hessAprox     = 'none', 'sr1', 'dfp', 'bfgs'
x0 = randn(2,1);
% x0 = [-3 10.234];
% x0 = [1.6534 0.74523];
trustRegion(f,g,x0,[],opt)

                            
                            
                            