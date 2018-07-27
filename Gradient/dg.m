function [ a ] = dg( z, fcn )
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
    for i = 1:length(z)
        if z(i) > 0
            a(i) = 0;
        else 
            a(i) = 1;
        end
    end
elseif strcmp(fcn,'relu2')
    a = max(0.01*z',z')';
    for i = 1:length(z)
        if z(i) > 0
            a(i) = 0.01;
        else 
            a(i) = 1;
        end
    end
elseif strcmp(fcn,'sigmoid')
    l = size(z,1);      %Number of lines
    c = size(z,2);      %Number of Columns

    %Calculate the sigmoid function for each term
    for i = 1:l
        for j = 1:c
            a(i,j) = exp(-z(i,j))/(1+exp(-z(i,j)))^2;
        end
    end
elseif strcmp(fcn,'softmax')
    disp('not implemented');
    return
else
    disp('Erro, o tipo da função não é reconhecido !');
end
    
end

