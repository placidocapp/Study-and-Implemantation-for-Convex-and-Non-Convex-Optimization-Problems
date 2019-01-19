clear all;
close all;
clc

%% Problem

%   min c'*x
%   subject to
%       A*x == b

c = [10; 8; 9; 5];
Aeq = [1 0 0 0
       0 1 0 0];
beq = [1; 0];
A = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];
b = [2;2;2;2];

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