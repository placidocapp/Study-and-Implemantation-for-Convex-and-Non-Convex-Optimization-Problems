function [ y ] = df( x )
%   Derivative of the function to solve
% There are some example functions to test, don't forget to change de
% function too ! (and the n value)

%% 1) Simple Quadratic  n = 1
% y = 2*x(1,1);

%% 2) Sin function  n = 1
% y = cos(x(1,1));

%% 3) Simple Quadratic  n = 1
% y = 2*x(1,1) + 3;

%% 4) Ackley N. 2 Function from http://benchmarkfcns.xyz/benchmarkfcns/ackleyn2fcn.html
   %n = 2       solutionn must be -200 and x=(0,0)
y(1,1) = (4*x(1).*exp(-(x(1).^2 + x(2).^2)^(1/2)/50))./(x(1).^2 + x(2).^2)^(1/2);
y(2,1) = (4*x(2).*exp(-(x(1).^2 + x(2).^2)^(1/2)/50))./(x(1).^2 + x(2).^2)^(1/2);
 
%% 5) Schwefel 2.20 Function from http://benchmarkfcns.xyz/benchmarkfcns/schwefel220fcn.html
   %n = 3       solutionn must be 0 and x=(0,0,0)
% y(1,1) = sign(x(1));
% y(2,1) = sign(x(2));
% y(3,1) = sign(x(3));
end

