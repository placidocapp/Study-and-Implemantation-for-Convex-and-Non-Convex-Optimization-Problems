
function [ J ] = fcost( w1, b1, w2, b2, w3, b3,...
    x, y, activation )
%% Initializations

m = size(x,2);    %Number of training examples
J = 0;
    
%% Calculate the cost function
for i = 1:m
    % Foward propagation
    z1 = w1*x(:,i) + b1;
    a1 = g(z1, activation);
    z2 = w2*a1 + b2;
    a2 = g(z2, activation);
    z3 = w3*a2 + b3;
    a3 = g(z3, 'sigmoid');
    
    %Cost function
    J = J - ( y(:,i).*log(a3) + (1-y(:,i)).*log(1-a3) );
end

J = sum(J)/m;




