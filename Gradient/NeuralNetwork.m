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
y_train = y_train(1:train_size);
y_test = y_test(1:test_size);
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

%%  Parameters

%For the main loop
ChooseMethod = 0;   %if 0 uses batch gradient, if 1 uses mini batch 
                    %gradient (note that if the mini batch has size 1 it's 
                    %a sthocastic gradient descent). if 2 use batch
                    %gradient with momentum.if 3 use desterov method
                    
maxIter = 1000;     %Stop creteria
grad_eps = 10^-10;  %Stop creteria
m = size(x_train,2);    %Number of training examples
n = size(x_train,1);    %Number of features
eps = 10^-8;

%The neural network have just 3 layers, but one can change the size of
%them
nlayers = 3;            %Number of layers
layers = [ n; 20; 20; 10 ]; %Size of each layer, the first layer is the                             
                            %number of features itself
activation = 'sigmoid';     %Activation function for all but the last layer 

%Gradients step length
alpha = 0.1;              %Step lenght

%Batch gradient [0, m]
batch_size = m;      

%In case we want batch gradient descent the batch_size is m
if ChooseMethod == 0
    batch_size = m;
end

%Momentum
lambda = 0.2;

%Adadelta
forget_rate = 0.8;

%RMSprop
beta = 0.5;

%Adam
beta1 = 0.9;
beta2 = 0.9;
                       
%% Initializations
rng(157024)
w1 = 0.1*randn(layers(2), layers(1));
w2 = 0.1*randn(layers(3), layers(2));
w3 = 0.1*randn(layers(4), layers(3));
b1 = 0.1*randn(layers(2),1);
b2 = 0.1*randn(layers(3),1);
b3 = 0.1*randn(layers(4),1);
J = zeros(maxIter,1);

%Variable for momentum
u1 = zeros(size(w1));
u2 = zeros(size(b1));
u3 = zeros(size(w2));
u4 = zeros(size(b2));
u5 = zeros(size(w3));
u6 = zeros(size(b3));

%Variables for best sol
best_J = 10^10;
best_w1 = zeros(size(w1));
best_w2 = zeros(size(w2));
best_w3 = zeros(size(w3));
best_b1 = zeros(size(b1));
best_b2 = zeros(size(b2));
best_b3 = zeros(size(b3));

%Variables for Adagrad
Gw1 = zeros(size(w1));
Gb1 = zeros(size(b1));
Gw2 = zeros(size(w2));
Gb2 = zeros(size(b2));
Gw3 = zeros(size(w3));
Gb3 = zeros(size(b3));

%Variables for Adadelta
Dw1 = rand(size(w1));
Db1 = rand(size(b1));
Dw2 = rand(size(w2));
Db2 = rand(size(b2));
Dw3 = rand(size(w3));
Db3 = rand(size(b3));
DDw1 = 0.1*rand(size(w1));
DDb1 = 0.1*rand(size(b1));
DDw2 = 0.1*rand(size(w2));
DDb2 = 0.1*rand(size(b2));
DDw3 = 0.1*rand(size(w3));
DDb3 = 0.1*rand(size(b3));

%Variables for Adam
mw1 = zeros(size(w1));
mb1 = zeros(size(b1));
mw2 = zeros(size(w2));
mb2 = zeros(size(b2));
mw3 = zeros(size(w3));
mb3 = zeros(size(b3));
vw1 = rand(size(w1));
vb1 = rand(size(b1));
vw2 = rand(size(w2));
vb2 = rand(size(b2));
vw3 = rand(size(w3));
vb3 = rand(size(b3));


%% Algorithm

tic
for k = 1:maxIter
    %% Batch sizes
    
    %Next Batch 
    batch = 1:m;
    batch = batch(randperm(m)); %Shuffle the batch
    batch = batch(1:batch_size);

    %% Function value and it's derivatives
    %Calculate the function value and derivatives with foward and back
    %propagation
    J(k) = fcost( w1, b1, w2, b2, w3, b3, x_train,...
        y_train, activation );
    
    %If the we need the nesterov method we need to choose the gradient for
    %a point ahead, otherwise we use regular gradient
    if ChooseMethod == 3
        [ dw1, db1, dw2, db2, dw3, db3 ] = fgrad( w1-lambda*u1...
            , b1-lambda*u2, w2-lambda*u3, b2-lambda*u4, w3-lambda*u5...
            , b3-lambda*u6, x_train,...
            y_train, activation, batch );
    else
        [ dw1, db1, dw2, db2, dw3, db3 ] = fgrad( w1, b1, w2, b2, w3, b3,...
            x_train ,y_train, activation, batch );
    end
    k
    kcost = J(k)      %Uncoment to see the evolution of cost

    %% Save the best solution
    if J(k) < best_J
        best_J = J(k);
        best_w1 = w1;
        best_w2 = w2;
        best_w3 = w3;
        best_b1 = b1;
        best_b2 = b2;
        best_b3 = b3;
    end
    
    %% Gradient methods
    if ChooseMethod == 0 || ChooseMethod == 1
        % Mini Batch Gratient descent
        w1 = w1 - alpha*dw1;
        b1 = b1 - alpha*db1;
        w2 = w2 - alpha*dw2;
        b2 = b2 - alpha*db2;
        w3 = w3 - alpha*dw3;
        b3 = b3 - alpha*db3;
        
        %One can decrease the alpha value on each iteration, it's good for
        %stochastic gradient descent
        cte = (alpha*0.1)^(-maxIter);
        alpha = alpha*cte;
    elseif ChooseMethod == 2
        % Momentum gradient descent
        %Update the momentum variable
        u1 = lambda*u1 + alpha*dw1;
        u2 = lambda*u2 + alpha*db1;
        u3 = lambda*u3 + alpha*dw2;
        u4 = lambda*u4 + alpha*db2;
        u5 = lambda*u5 + alpha*dw3;
        u6 = lambda*u6 + alpha*db3;
        
        %Update the weights points
        w1 = w1 - u1;
        b1 = b1 - u2;
        w2 = w2 - u3;
        b2 = b2 - u4;
        w3 = w3 - u5;
        b3 = b3 - u6;
    elseif ChooseMethod == 3
        %Nesterov Gradient descent
        %Update the momentum variable
        u1 = lambda*u1 + alpha*dw1;
        u2 = lambda*u2 + alpha*db1;
        u3 = lambda*u3 + alpha*dw2;
        u4 = lambda*u4 + alpha*db2;
        u5 = lambda*u5 + alpha*dw3;
        u6 = lambda*u6 + alpha*db3;
        
        %Update
        w1 = w1 - u1;
        b1 = b1 - u2;
        w2 = w2 - u3;
        b2 = b2 - u4;
        w3 = w3 - u5;
        b3 = b3 - u6;
    elseif ChooseMethod == 4
        %Adagrad Gradient Descent
        %Auxiliar variable G
        Gw1 = Gw1*forget_rate + dw1.^2;
        Gb1 = Gb1*forget_rate + db1.^2;
        Gw2 = Gw2*forget_rate + dw2.^2;
        Gb2 = Gb2*forget_rate + db2.^2;
        Gw3 = Gw3*forget_rate + dw3.^2;
        Gb3 = Gb3*forget_rate + db3.^2;
        
        %Update
        w1 = w1 - alpha./sqrt(Gw1+eps*ones(size(dw1))).*dw1;
        b1 = b1 - alpha./sqrt(Gb1+eps*ones(size(db1))).*db1;
        w2 = w2 - alpha./sqrt(Gw2+eps*ones(size(dw2))).*dw2;
        b2 = b2 - alpha./sqrt(Gb2+eps*ones(size(db2))).*db2;
        w3 = w3 - alpha./sqrt(Gw3+eps*ones(size(dw3))).*dw3;
        b3 = b3 - alpha./sqrt(Gb3+eps*ones(size(db3))).*db3;
    elseif ChooseMethod == 5
        %Adadelta Gradient Descent
        %Auxiliar variable G
        Gw1 = Gw1*forget_rate + (1-forget_rate)*dw1.^2;
        Gb1 = Gb1*forget_rate + (1-forget_rate)*db1.^2;
        Gw2 = Gw2*forget_rate + (1-forget_rate)*dw2.^2;
        Gb2 = Gb2*forget_rate + (1-forget_rate)*db2.^2;
        Gw3 = Gw3*forget_rate + (1-forget_rate)*dw3.^2;
        Gb3 = Gb3*forget_rate + (1-forget_rate)*db3.^2;
        
        %Update Delta
        Dw1 = - ((DDw1+eps*ones(size(dw1))).^0.5)./((Gw1+eps*ones(size(dw1))).^0.5).*dw1;
        Db1 = - ((DDb1+eps*ones(size(db1))).^0.5)./((Gb1+eps*ones(size(db1))).^0.5).*db1;
        Dw2 = - ((DDw2+eps*ones(size(dw2))).^0.5)./((Gw2+eps*ones(size(dw2))).^0.5).*dw2;
        Db2 = - ((DDb2+eps*ones(size(db2))).^0.5)./((Gb2+eps*ones(size(db2))).^0.5).*db2;
        Dw3 = - ((DDw3+eps*ones(size(dw3))).^0.5)./((Gw3+eps*ones(size(dw3))).^0.5).*dw3;
        Db3 = - ((DDb3+eps*ones(size(db3))).^0.5)./((Gb3+eps*ones(size(db3))).^0.5).*db3;
        
        %Auxiliar variable D
        DDw1 = DDw1*forget_rate + (1-forget_rate)*Dw1.^2;
        DDb1 = DDb1*forget_rate + (1-forget_rate)*Db1.^2;
        DDw2 = DDw2*forget_rate + (1-forget_rate)*Dw2.^2;
        DDb2 = DDb2*forget_rate + (1-forget_rate)*Db2.^2;
        DDw3 = DDw3*forget_rate + (1-forget_rate)*Dw3.^2;
        DDb3 = DDb3*forget_rate + (1-forget_rate)*Db3.^2;
        
        %Final variable
        w1 = w1 + Dw1;
        b1 = b1 + Db1;
        w2 = w2 + Dw2;
        b2 = b2 + Db2;
        w3 = w3 + Dw3;
        b3 = b3 + Db3;
        
    elseif ChooseMethod == 6
        %RMSprop Gradient descent
        %Update the RMS
        u1 = beta*u1 + (1-beta)*dw1.^2;
        u2 = beta*u2 + (1-beta)*db1.^2;
        u3 = beta*u3 + (1-beta)*dw2.^2;
        u4 = beta*u4 + (1-beta)*db2.^2;
        u5 = beta*u5 + (1-beta)*dw3.^2;
        u6 = beta*u6 + (1-beta)*db3.^2;
        
        %Update 
        w1 = w1 - alpha*dw1./sqrt(u1+eps*ones(size(dw1)));
        b1 = b1 - alpha*db1./sqrt(u2+eps*ones(size(db1)));
        w2 = w2 - alpha*dw2./sqrt(u3+eps*ones(size(dw2)));
        b2 = b2 - alpha*db2./sqrt(u4+eps*ones(size(db2)));
        w3 = w3 - alpha*dw3./sqrt(u5+eps*ones(size(dw3)));
        b3 = b3 - alpha*db3./sqrt(u6+eps*ones(size(db3)));
    elseif ChooseMethod == 7
        %Adam Gradient Descent
        %First Moment update
        mw1 = beta1*mw1 + (1-beta1)*dw1;
        mb1 = beta1*mb1 + (1-beta1)*db1;
        mw2 = beta1*mw2 + (1-beta1)*dw2;
        mb2 = beta1*mb2 + (1-beta1)*db2;
        mw3 = beta1*mw3 + (1-beta1)*dw3;
        mb3 = beta1*mb3 + (1-beta1)*db3;
        
        %Second Moment update
        vw1 = beta2*vw1 + (1-beta2)*dw1.^2;
        vb1 = beta2*vb1 + (1-beta2)*db1.^2;
        vw2 = beta2*vw2 + (1-beta2)*dw2.^2;
        vb2 = beta2*vb2 + (1-beta2)*db2.^2;
        vw3 = beta2*vw3 + (1-beta2)*dw3.^2;
        vb3 = beta2*vb3 + (1-beta2)*db3.^2;
        
        %Update
        w1 = w1 - alpha*(mw1/(1-beta1^k))./sqrt(vw1/(1-beta2^k)+eps*ones(size(dw1)));
        b1 = b1 - alpha*(mb1/(1-beta1^k))./sqrt(vb1/(1-beta2^k)+eps*ones(size(db1)));
        w2 = w2 - alpha*(mw2/(1-beta1^k))./sqrt(vw2/(1-beta2^k)+eps*ones(size(dw2)));
        b2 = b2 - alpha*(mb2/(1-beta1^k))./sqrt(vb2/(1-beta2^k)+eps*ones(size(db2)));
        w3 = w3 - alpha*(mw3/(1-beta1^k))./sqrt(vw3/(1-beta2^k)+eps*ones(size(dw3)));
        b3 = b3 - alpha*(mb3/(1-beta1^k))./sqrt(vb3/(1-beta2^k)+eps*ones(size(db3)));
    elseif ChooseMethod == 8
        %Adamax Gradient Descent
        %First Moment update
        mw1 = beta1*mw1 + (1-beta1)*dw1;
        mb1 = beta1*mb1 + (1-beta1)*db1;
        mw2 = beta1*mw2 + (1-beta1)*dw2;
        mb2 = beta1*mb2 + (1-beta1)*db2;
        mw3 = beta1*mw3 + (1-beta1)*dw3;
        mb3 = beta1*mb3 + (1-beta1)*db3;
        
        %Second Moment update
        vw1 = max(beta2*vw1, abs(dw1));
        vb1 = max(beta2*vb1, abs(db1));
        vw2 = max(beta2*vw2, abs(dw2));
        vb2 = max(beta2*vb2, abs(db2));
        vw3 = max(beta2*vw3, abs(dw3));
        vb3 = max(beta2*vb3, abs(db3));
        
        %Update
        w1 = w1 - alpha*(mw1/(1-beta1^k))./(vw1+eps*ones(size(dw1)));
        b1 = b1 - alpha*(mb1/(1-beta1^k))./(vb1+eps*ones(size(db1)));
        w2 = w2 - alpha*(mw2/(1-beta1^k))./(vw2+eps*ones(size(dw2)));
        b2 = b2 - alpha*(mb2/(1-beta1^k))./(vb2+eps*ones(size(db2)));
        w3 = w3 - alpha*(mw3/(1-beta1^k))./(vw3+eps*ones(size(dw3)));
        b3 = b3 - alpha*(mb3/(1-beta1^k))./(vb3+eps*ones(size(db3)));
        
    end
    
    %% Stop creteria
%     norm([norm(dw1) norm(db1) norm(dw2) norm(db2) norm(dw3) norm(db3)])
    if norm([norm(dw1) norm(db1) norm(dw2) norm(db2) norm(dw3) norm(db3)])...
            <= grad_eps
        break;
    end
    
    %% Plot 
    figure(1), plot(J(1:k),'b'), hold on
end
time = toc

%Before prediction uses the best values achived
w1 = best_w1;
w2 = best_w2;
w3 = best_w3;
b1 = best_b1;
b2 = best_b2;
b3 = best_b3;

%% Prediction

y_pred = zeros(test_size,1);

for i = 1:test_size
    
    %Foward propagation
    z1 = w1*x_test(:,i) + b1;
    a1 = g(z1, activation);
    z2 = w2*a1 + b2;
    a2 = g(z2, activation);
    z3 = w3*a2 + b3;
    a3 = g(z3, 'sigmoid');

    %Transform in numbers again
    [aux, val] = max(a3);
    y_pred(i) = val - 1;

end

correct = sum(y_pred == y_test);
accuracy = correct/test_size

save('Test','J','time','accuracy','best_J')