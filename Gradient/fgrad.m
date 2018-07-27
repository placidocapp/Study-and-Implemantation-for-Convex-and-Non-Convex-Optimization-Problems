function [ dw1, db1, dw2, db2, dw3, db3 ] = fgrad( w1, b1, w2, b2, w3, b3,...
    x, y, activation, batch )
%% Initializations

batchSize = length(batch);
dw1 = zeros();
db1 = zeros();
dw2 = zeros();
db2 = zeros();
dw3 = zeros();
db3 = zeros();

%Use all examples available is the default
if ~exist('batch','var')
    batch = 1:m;
end

%% Algorithm

for k = 1:batchSize
    %Go to the indices sorted in batch
    i = batch(k);
    
    %% Foward propagation
    
    z1 = w1*x(:,i) + b1;
    a1 = g(z1, activation);
    z2 = w2*a1 + b2;
    a2 = g(z2, activation);
    z3 = w3*a2 + b3;
    a3 = g(z3, 'sigmoid');

    %% Back Propagation

    %Last layer
%     da3 = - y(:,i)./a3 + (1-y(:,i))./(1-a3);
%     dz3 = da3.*dg(z3,'sigmoid');
    %We can simplify the 2 steps above multiplying them and we obtain the
    %below equation
    dz3 = a3 - y(:,i);
    dw3 = dw3 + dz3*a2';
    db3 = db3 + dz3;
    
    %Second layer
    da2 = w3'*dz3;
    dz2 = da2.*dg(z2,activation);
    dw2 = dw2 + dz2*a1';
    db2 = db2 + dz2;
    
    %First layer
    da1 = w2'*dz2;
    dz1 = da1.*dg(z1,activation);
    dw1 = dw1 + dz1*x(:,i)';
    db1 = db1 + dz1;
    
end

dw1 = dw1/batchSize;
db1 = db1/batchSize;
dw2 = dw2/batchSize;
db2 = db2/batchSize;
dw3 = dw3/batchSize;
db3 = db3/batchSize;

end

