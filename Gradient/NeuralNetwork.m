close all;
clear all;
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

%The logistic regression can only classify into two different groups, so i
%will modify the original y 
train_size = 1000;
test_size = 1000;
y_train = (y_train(1:train_size) == 4);
y_test = (y_test(1:test_size) == 4);
x_train = x_train(:,1:train_size);
x_test = x_test(:,1:test_size);

%% Change the format of y

%Here we need that for each example the y variable become a column with
%10 values, just one for the respective number (from 0 to 9) and zeros
%for the others
aux = zeros(10,train_size);
for i = 1:train_size
    aux(:,i) = (0:9)'==y_train(i);
end
y_train = aux;

aux2 = zeros(10,test_size);
for i = 1:test_size
    aux(:,i) = (0:9)'==y_test(i);
end
y_test = aux;

%%  Parameters

maxIter = 1000;
m = size(x_train,2);    %Number of training examples
n = size(x_train,1);    %Number of features
alpha = 0.1;              %Step lenght
nlayers = 3;            %Number of layers
layers = [ n; 20; 10 ]; %Size of each layer, the first layer is the                             
                        %number of features itself
activation = 'sigmoid';    %Activation function for all but the last layer 
                        
%% Initializations

w1 = 0.1*randn(layers(2), layers(1));
w2 = 0.1*randn(layers(3), layers(2));
b1 = 0.1*randn(layers(2),1);
b2 = 0.1*randn(layers(3),1);

%% Algorithm

for k = 1:maxIter
    
    %Calculate the function value and derivatives with foward and back
    %propagation
    J = 0;
    for i = 1:m
        [a2, dw1, db1, dw2, db2] = fpropag( w1, b1, w2, b2, x_train(:,i),...
            y_train(:,i), activation );
        J = J - ( y_train(i)*log(a2) + (1-y_train(i))*log(1-a2) );
    end
    J = sum(J)/m
    
    %Gratient descent
%     w1 = w1 - alpha*dw1;
%     b1 = b1 - alpha*db1;
    w2 = w2 - alpha*dw2;
    b2 = b2 - alpha*db2;
    
end












