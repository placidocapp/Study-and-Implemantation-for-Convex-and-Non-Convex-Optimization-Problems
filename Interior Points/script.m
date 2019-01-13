clear all;
close all;
clc

%% Problem

%   min c'*x
%   subject to
%       A*x == b

c = [10; 8; 9; 5];
A = [0 2 1 0
     0 2 1 6
     1 0 0 1];
b = [4; 8; 6];

%Solution
cvx_begin quiet
variables x(4)
minimize(c'*x)
subject to
    A*x == b
    x >= 0
cvx_end

cvx_status
x
cvx_optval

%% Call lp

[x_opt,f_opt,status] = lp(c,A,b,q,1)