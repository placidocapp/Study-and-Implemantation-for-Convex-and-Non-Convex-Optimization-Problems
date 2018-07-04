% Algorithm from Numerical Optimization. Jorge Nocedal and Stephen J. Wright.
% Second edition. Zoom algorithm can be found at pages 69.
% MATLAB code by Pl�cido Campos
% Last modified: july 4, 2018

clear all;
close all;
clc

%%  Parameters

maxIter = 100;  %Maximun number of iterations
n = 2;          %Dimension of the problem   
eta = 0.20;     %Parameter that decides if the algorithm steps or not,
                %the book recomends it between (0,0.25)

%% Inicialization

p = zeros(maxIter,n);       %Direction of the step
ro = zeros(maxIter,1);      %rate of real decrease divided by estimatede decrease
x = zeros(maxIter,n);       %Actual point
x(1,:) = [10.3 17.2];
delta = zeros(maxIter,1);   %Radius of trust area
delta(1) = 1;
deltaMax = 10;

%% Function

%Valores da fun��o
A = randn(2);
A = A'*A;
A = A'+ A;
b = randn(2,1);

f = @(x) 0.5*x'*A*x + b'*x;
g = @(x) A*x + b;
H = @(x) A;

m = @(x,p) f(x) + g(x)'*p + 0.5*p'*B*p;

%% Algorithm

for k = 1:maxIter
    %% Solve argmin mk(p) to find pk
    %Se for semidefinido positivo
    aux = g(x)'*B(:,:,k)*g(x);
    if aux >= 0
        tau(k) = 1;
    else 
        tau(k) = min(1, norm(g(x))^3/(delta(k)*aux) );
    end
    %% Decide to step or to rise/decrease the area delta
    
    %Evaluate ro(k)
    ro(k) = ( f(x(k)) - f(x(k)+p(k)) )/( m(0) - m(p(k)) );
    if ro(k) < 0.25
        delta(k+1) = 0.25*delta(k);
    else
        if (ro(k) > 3/4 || norm(p(k)) == delta)
            delta(k+1) = min(2*delta(k), deltaMax);
        else
            delta(k+1) = delta(k);
        end
    end
    
    %Decide to take the step or not
    if ro(k) > eta
        x(k+1) = x(k) + p(k);
    else
        x(k+1) = x(k);
    end
end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        