function [ a2, dw1, db1, dw2, db2 ] = fpropag( w1, b1, w2, b2, x, y, activation )
%Foward Propagation

%% Initializations
% z1 = zeros(layers(2),1);
% z2 = zeros(layers(3),1);
% a1 = zeros(layers(2),1);
% a2 = zeros(layers(3),1);
m = size(x,2);    %Number of training examples

%% Foward propagation

z1 = w1*x + b1;
a1 = g(z1, activation);
z2 = w2*a1 + b2;
a2 = g(z2, 'sigmoid');

%% Back Propagation

%Last layer
da2 = - y./a2 + (1-y)./(1-a2);
dz2 = da2.*dg(z2,'sigmoid');
% dz2 = a2 - y;
dw2 = dz2*a1';
db2 = dz2;

%First layer
da1 = w2'*dz2;
dz1 = da1.*dg(z1,activation);
dw1 = dz1*x';
db1 = dz1;


