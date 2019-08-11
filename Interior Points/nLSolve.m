function [x_opt,f_opt,status] = nLSolve(f,fx,fxx,h,hx,hxx,g,gx,gxx,x0,opt)
%This function aims to solve an optimization problem of the form
%       min f(x)
%       subject to 
%           h(x) =  0
%           g(x) >= 0
%The arguments of this function are:
%       f:   objective function
%       fx:  gradient of f
%       fxx: Hessian of f (optional)
%       h:   Equality constraints
%       hx:  gradient of h
%       hxx: Hessian of h (optional)
%       g:   Inequality constraints
%       gx:  gradient of f
%       gxx: Hessian of f (optional)
%       x0:  initial condition
%       opt: options for algorithm to be used

%% Arguments

hessianAvailable = 1;
if ~exist('fxx','var')
    hessianAvailable = 0;
end

n = size(x0,1)

%% Initial Solution and Inicialization

s0 = rand(n,1)

end

