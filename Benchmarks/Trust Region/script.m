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
maxIter = 1000;

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
    opt = options('dogleg','sr1');
    %algorithm     = 'cauchy', 'dogleg', 'subproblem'
    %hessAprox     = 'none', 'sr1', 'dfp', 'bfgs'
    rng(1);
    for j = 1:10
        x0 = 10*randn(2,1);
        [a, b, c, d, e] = trustRegion(f,g,x0,B,opt);
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

%% Chama algoritmo


[f,g,B,sol,fsol] = benchmark(10);
opt = options('cauchy','bfgs');
%algorithm     = 'cauchy', 'dogleg', 'subproblem'
%hessAprox     = 'none', 'sr1', 'dfp', 'bfgs'
x0 = randn(2,1);
% x0 = [-3 10.234];
% x0 = [1.6534 0.74523];
trustRegion(f,g,x0,B,opt)

                            
                            
                            