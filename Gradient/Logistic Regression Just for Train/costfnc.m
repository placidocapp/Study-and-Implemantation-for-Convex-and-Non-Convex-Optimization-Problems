function [ J, dw, db ] = costfnc( w, b, x, y )
%Calculate the cost function for the logistic regression and it's
%derivatives

%% Parameters

m = size(x,2);    %Number of training examples
n = size(x,1);    %Number of features

%% Inicializations

J = zeros(m,1);   %Cost function
z = zeros(m,1);        %Linear regression for y
dw = zeros(n,1);        %Derivative of J with respect to w
db = zeros(1);        %Derivative of J with respect to b
dz = zeros(m,1);          %derivative of sigmoi with respect to z
a = zeros(m,1);         %Sigmoid value

for i = 1:m
    %For all features we calculate each example linear's regression
    z(i) = w'*x(:,i) + b;
    
    %For all features we calculate each example sigmoid's function
    a(i) = sigmoid(z(i));
    
    %Increment the cost function value
    J(i) = J(i) - ( y(i)*log(a(i)) + ( 1-y(i) )*log( 1-a(i) ) );
    
    %Calculate the derivatives for this example
    dz(i) = a(i) - y(i);
    dw(:,1) = dw(:,1) + x(:,i)*( dz(i) ); 
    db = db + dz(i);
end

J = sum(J)/m;
dw = dw/m;
db = db/m;

end

