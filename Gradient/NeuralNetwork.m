close all;
cleal all;
clc

%% Load dataset

disp('Loading data ...')

%We recive here the pictures in x and labels in y. The variables are
%x_train, y_train, x_test and y_test. For trainning
load('MNISTdataset.mat')

% Now to work on the data I will transform each pixel in one feature, so I
% will change the format of the matrices x_train and x_test
x_train = reshape(x_train,[400 60000]);
x_test = reshape(x_test,[400 10000]);

%%  Parameters

maxIter = 1000;
m = size(x_train,2);    %Number of training examples
n = size(x_train,1);    %Number of features
alpha = 1;              %Step lenght
nlayers = 3;            %Number of layers
layers = [ n; 20; 10 ]; %Size of each layer, the first layer is the number
                        %of features itself

%% Initializations

%Weights for activation functions
w1 = randn(layers(2), layers(1));
w2 = randn(layers(3), layers(2));
b1 = randn(layers(2),1);
b2 = randn(layers(3),1);
z1 = zeros(layers(2),1);
z2 = zeros(layers(3),1);
a1 = zeros(layers(2),1);
a2 = zeros(layers(3),1);

%% Algorithm

%Foward propagation
















