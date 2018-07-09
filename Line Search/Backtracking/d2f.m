function [ y ] = d2f( x )
%   Derivative of the function to solve
% There are some example functions to test, don't forget to change de
% function too ! (and the n value)

%% 1) Simple Quadratic  n = 1
% y = 2;

%% 2) Sin function  n = 1
y = -sin(x(1,1));

%% 3) Simple Quadratic  n = 1
% y = 2;

%% 4) Ackley N. 2 Function from http://benchmarkfcns.xyz/benchmarkfcns/ackleyn2fcn.html
   %n = 2       solutionn must be -200 and x=(0,0)
y = [ (4*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(1/2) - (2*x(1)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2),                                                    - (2*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2);
                                                    - (2*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2), (4*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(1/2) - (2*x(2)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(2)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2)] ;
%% 5) Schwefel 2.20 Function from http://benchmarkfcns.xyz/benchmarkfcns/schwefel220fcn.html
   %n = 3       solutionn must be 0 and x=(0,0,0)
% y(1,1) = 0;
% y(2,1) = 0;
% y(3,1) = 0;

end

