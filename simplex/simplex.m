function [x_opt,f_opt] = simplex(c,A,b,x0,base)
%simpex is na algorithm for solving linear problems
%Here A is already made of equality constrains
%       min c'*x
%       subject to
%           A*x = b

%% Variables

m = size(A,1);
n = size(A,2);

if m ~= length(b)
    disp('Number of elements in b and lines in A must be equal');
    return
end

if n ~= length(c)
    A
    c
    disp('Number of elements in c and columns in A must be equal');
    return
end

tableau = [c'*x0,    c';
           x0(base), A]

%% ALgorithm

while(1)
    

x_opt = [];
f_opt = [];
end