clear all;
close all;
clc

%% Problem

%   min c'*x
%   subject to
%       A*x   <= b
%       Aeq*x == beq

c = -[10; 12; 12];
A = [1 2 2;
     2 1 2;
     2 2 1];
b = [20; 20; 20];
Aeq = [];
beq = [];

%Solution
cvx_begin quiet
variables x(3)
minimize(c'*x)
subject to
    A*x <= b
%     Aeq*x == beq
    x >= 0
cvx_end

cvx_status
x
cvx_optval

%% Call lp

[x_opt,f_opt] = lp(c,A,b,Aeq,beq)