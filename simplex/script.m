clear all;
close all;
clc

%% Problem

%   min c'*x
%   subject to
%       A*x   <= b
%       Aeq*x == beq

c = -[1; 1; 1];
Aeq = [1 2 3;
     -1 2 6;
     0 4 9];
beq = [3; 2; 5];
A = [0 0 3];
b = [1];

%Solution
cvx_begin quiet
variables x(3)
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