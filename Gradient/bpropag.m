function [ dw1, db1, dw2, db2 ] = ...
    bpropag( w1, b1, w2, b2, x, y, a2, activation )
%Do backpropagation to find the derivatives of J with respect to the
%wieghts in the activation functions

%Derivative of J with respect to a2 (for just one example) 
da2 = - y./a2 + (1-y)./(1-a2);
dw2 = x*dg(z2)

end

