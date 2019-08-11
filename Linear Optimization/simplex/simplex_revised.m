function [x_opt,f_opt,status] = simplex_revised(c,A,b,Aeq,beq);
%Revised simplex is na algorithm for solving linear problems
%Here A is already made of equality constrains
%       min c'*x
%       subject to
%           A*x <= b
%           Aeq*x = beq

%% Variables

m = size(A,1) + size(Aeq,1);    %number of constrains
n = size(A,2);                  %number of variables
na = size(A,1);                 %number of slack variables
ny = m - size(A,1);             %number of extra variables

%Get together all restrictions
A = [A; Aeq];
b = [b; beq];

if m ~= length(b)
    disp('Number of elements in b and lines in A must be equal');
    return
end

if n ~= length(c)
    disp('Number of elements in c and columns in A must be equal');
    return
end

%% Fist Phase - Find an Feasible Solution

%If there are not enought slack variables, than use the extra variables 
%to find an initial solution
if ny > 0
    %Slack variables and aditional variables
    A = [A, eye(m)];
    c = [c; zeros(na,1)];
    n = n + ny + na;

    %Initial solution to auxiliary problem
    x0 = [zeros(n,1); b];

    %Vector with base variables
    base = (n+1):(n+m);
    B = eye(m);
    
    %Auxiliary c for extra variables
    caux = c;
    c = [zeros(n-ny+m,1); ones(ny,1)];

    %Call the solver
    simplex_revised();
    
end
    
end

