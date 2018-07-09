function [ y ] = f( x )
%   Function to solve
% There are some example functions to test, don't forget to change de
% derivative too ! (and the n value)

%% 1) Simple Quadratic  n = 1
% y = x(1,1)^2;

%% 2) Sin function  n = 1
% y = sin(x(1,1));

%% 3) Simple Quadratic  n = 1
% y = x(1,1)^2 + 3*x(1,1) + 50;

%% 4) Ackley N. 2 Function from http://benchmarkfcns.xyz/benchmarkfcns/ackleyn2fcn.html
   %n = 2       solutionn must be -200 and x=(0,0)
y = -200 * exp(-0.02 * sqrt((x(1,1) .^ 2) + (x(2,1) .^ 2)));

%% 5) Schwefel 2.20 Function from http://benchmarkfcns.xyz/benchmarkfcns/schwefel220fcn.html
   %n = 3       solutionn must be 0 and x=(0,0,0)
% y = sum(abs(x), 1);

end

