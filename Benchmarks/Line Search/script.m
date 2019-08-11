% MATLAB code by Plácido Campos based on Jorge Nocedal and Stephen J. Wright.
% Numerical Optimization, Second Edition.

clear all;
close all;
clc

format long;
%% Parameters 

% Wolfe conditions constants
c1 = 10^-4;        
c2 = 0.5;

%Stop criteria for norm of gradient
eps = 10^-8;

%Backtracking division constant
ro = 0.5;

%Limit of iterations
maxIter = 100;

%Size of the problem (Selects size based on the gradient size)
n = 2;

kfinal = -1;        %The iteration that the algorithm stopped
stepMethod = 1;     %if 0 uses backtracking to choose the step lenght else 
                    %uses the zoom algorithm

%% Function
% fopt = zeros(5, 2);
% xopt = zeros(5, 2, 2);
% met = {'newton'; 'bfgs'; 'dfp'; 'sr1'; 'gradient'};
% 
% for i = 1:10
%     [f,g,B,sol] = benchmark(i);
%     for l = 1:5
%         k = 1;
%         for alpha = 'backtacking'
%             for modHess = 'cholesky'
%                 x0 = [10.13455 -5.3234];
%                 opt = options(string(met(l)),'backtacking','cholesky');
%                 [xopt(l,:,i), fopt(l,i)] = lineSearch(f,g,x0,B,opt,10^-8,1000)
%             end
%             k = k + 1;
%         end
%     end
% end

%% Chama algoritmo
xopt = zeros(10,1);
fopt = zeros(10,1);
time = zeros(10,1);
gnorm = zeros(10,1);
iter = zeros(10,1);
errox = zeros(10,1);
errof = zeros(10,1);

for i = 1:10
    [f,g,B,sol,fsol] = benchmark(i);
    opt = options('gradient','wolfeCond','cholesky');
    %algorithm     = 'newton', 'bfgs', 'dfp', 'sr1', 'gradient'
    %stepMethod    = 'backtacking','wolfeCond','bisection'
    %modHessMethod = 'norm2', 'cholesky'
    rng(1);
    for j = 1:10
        x0 = 10*randn(2,1);
%     x0 = [2.269 -1.549];
        [a, b, c, d, e] = lineSearch(f,g,x0,B,opt,10^-8,100);
        xopt(i) = xopt(i) + norm(a-sol)/10;
        fopt(i) = fopt(i) + abs(b-fsol)/10;
        time(i) = time(i) + c/10;
        gnorm(i) = gnorm(i) + d/10;
        iter(i) = iter(i) + e/10;
    end
    errox(i) = xopt(i);
    errof(i) = fopt(i);
end

T = table(errox,errof,time,gnorm,iter)
writetable(T,'data.xlsx','Sheet',1,'Range','A1')
 
return
%% Teste simples            

[f,g,B,sol,fsol] = benchmark(3);
opt = options('bfgs','wolfeCond','cholesky');
[xopt, fopt] = lineSearch(f,g,x0,B,opt,10^-8,1000)

