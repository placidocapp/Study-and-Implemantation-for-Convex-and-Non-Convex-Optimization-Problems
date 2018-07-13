function [ a ] = g( z, fcn )
%Calculate an activation function
%fcn can be 'relu', 'relu2' 'sigmoid'

%Relu is the default option
if ~exist('fcn','var')
    fcn = 'relu';
end

%Inicialization
a = zeros(size(z));

%Function calculation
if strcmp(fcn,'relu')
    a = max(0,z')';
elseif strcmp(fcn,'relu2')
    a = max(0.01*z',z')';
elseif strcmp(fcn,'sigmoid')
    l = size(z,1);      %Number of lines
    c = size(z,2);      %Number of Columns

    %Calculate the sigmoid function for each term
    for i = 1:l
        for j = 1:c
            a(i,j) = 1/(1+exp(-z(i,j)));
        end
    end
elseif strcmp(fcn,'softmax')
    l = size(z,1);      %Number of lines
    c = size(z,2);      %Number of Columns
    
    %Calculate the softmax function for each term
    for i = 1:l
        for j = 1:c
            a(i,j) = exp(z(i,j))/sum(exp(z(:,j)));
        end
    end
else
    disp('Erro, o tipo da função não é reconhecido !');
end
    
end

