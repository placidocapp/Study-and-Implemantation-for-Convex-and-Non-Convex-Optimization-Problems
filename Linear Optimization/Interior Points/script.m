clear all;
close all;
clc

%% Problem

%   min f(x)
%   subject to
%       h(x) = 0
%       g(x) >= 0

g = @(x) x^2;
h = @(x) 2*x + 3;



return
%Solution
cvx_begin quiet
variables x(4)
minimize(c'*x)
subject to
    A*x <= b
    Aeq*x == beq
    x >= 0
cvx_end

cvx_status
x
cvx_optval

%% Call lp

[x_opt,f_opt,status] = lp(c,A,b,Aeq,beq)