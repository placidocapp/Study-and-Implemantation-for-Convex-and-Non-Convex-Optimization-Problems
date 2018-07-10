% Algorithm from Numerical Optimization. Jorge Nocedal and Stephen J. Wright.
% Second edition. Zoom algorithm can be found at pages 60 and 61.
% MATLAB code by Plácido Campos
% Last modified: july 4, 2018

clear all;
close all;
clc

format long;
%% Parameters 

% Wolfe conditions
c1 = 10^-4;
c2 = 0.5;
alphaMax = 10;
eps = 10^-8;

%Backtracking
ro = 0.5;

%Todos
maxIter = 100;

%Size (Selects size based on the gradient size)
n = 2;

kfinal = -1;        %The iteration that the algorithm stopped
flagChooseStep = 1; %if 0 uses backtracking to choose the step lenght else 
                    %uses the zoom algorithm

%% Function

%Function values
A = randn(n);
A = A'*A;
A = A'+ A;
b = randn(n,1);

%Calculate the function value and its derivatives
% f = @(x) 0.5*x'*A*x + b'*x;
% g = @(x) A*x + b;
% B = @(x) A;
% sol = - (pinv(A)*b)';

f = @(x) sin(x(1))+sin(x(2));
g = @(x) [cos(x(1));
          cos(x(2))];
B = @(x) [  -sin(x(1)) 0;
            0          -sin(x(2))];
sol = [-1.570796327268957  -1.570796324699814];

% f = @(x) -200*exp(-0.2*sqrt(x(1)^2+x(2)^2));
% g = @(x) [(4*x(1).*exp(-(x(1).^2 + x(2).^2)^(1/2)/50))./(x(1).^2 + x(2).^2)^(1/2);
%           (4*x(2).*exp(-(x(1).^2 + x(2).^2)^(1/2)/50))./(x(1).^2 + x(2).^2)^(1/2)];
% B = @(x) [ (4*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(1/2) - (2*x(1)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2),                                                    - (2*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2);
%                                                     - (2*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2), (4*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(1/2) - (2*x(2)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(2)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2)] ;
% sol = [0 0];

%Rosenbrock
% a = 1;
% b = 100;
% f = @(x) (a-x(1))^2 + b*(x(2)-x(1)^2)^2; 
% g = @(x) [-2*(a - x(1))-4*b*x(1)*(x(2)-x(1)^2);
%           2*b*(x(2)-x(1)^2)];
% B = @(x) [2-4*b*x(2)+12*b*x(1)^2 -4*b*x(1);
%           -4*b*x(1)                2*b];       
% sol = [a a^2];


%% Backtraking Algorithm (Modified Hessian with 2-norm)

%Inicialization of x
x = zeros(maxIter,n);
x(1,:) = randn(1,n);

for flag = 0:2
    %Backtraking loop
    for k = 1:maxIter

        %Calculates the hessian, if some eigenvalue of H is negative than
        %change the hessian for itself plus eye*|most negative eigenvalue + lambda|
        H = B(x(k,:)');
        if flag == 0
            if sum(eig(H) < 0) > 0
                lambda = min(eig(H));
                H = H + eye*(-lambda+0.1);
            end
            d = -inv(H)*g(x(k,:)');
        elseif flag == 1
            if sum(eig(H) < 0) > 0
                [L,D] = mcfac(H);
                H = L*D*L';
            end
            d = -inv(H)*g(x(k,:)');
        else 
            d = -g(x(k,:)');
        end

        %Choose the step lenght
        if flagChooseStep == 0
            alpha = backtrack( f, g, x(k,:)', d, c1, ro );
        else
            alpha = lineSearch( f, g, x(k,:)', d, c1, c2 );
        end
        
        % Stop criteria
        if norm(g(x(k,:)')) < eps
            disp('Stopped at')
            kfinal = k
            break
        end

        %Find the next point
        x(k+1,:) = (x(k,:)' + alpha*d)';
    end

    if kfinal == -1
        kfinal = maxIter;
    end
    
    x(kfinal,:)
    norm(g(x(kfinal,:)'))
    f(x(kfinal,:)')

    %% Plot some graphs

    x_error = zeros(kfinal,1);
    f_error = zeros(kfinal,1);
    g_evol = zeros(kfinal,1);
    for i = 1:kfinal
        x_error(i) = norm(sol) - norm(x(i,:));
        f_error(i) = f(sol') - f(x(i,:)');
        g_evol(i) = norm(g(x(i,:)'));
    end

    k = 1:kfinal;

    subplot(3,1,1), plot(k,x_error), hold on, title('Error in position x')
    legend('Modified B with 2-norm','Modified B sith cholesky','Steepest descent')
    subplot(3,1,2), plot(k,f_error), hold on, title('Error in function value')
    legend('Modified B with 2-norm','Modified B sith cholesky','Steepest descent')
    subplot(3,1,3), plot(k,g_evol), hold on, title('Gradient evolution')
    legend('Modified B with 2-norm','Modified B sith cholesky','Steepest descent')
end
return
%% Backtraking Algorithm (Modified Hessian with cholesky factorization)

%Use this to find de derivative of a new function
% syms x;
% f = f(x);
% df = diff(f);

%Inicialization of x
x = zeros(maxIter,n);
x(1,:) = randn(1,n);
% x(1,:) = 0.333510833065806;

%Backtraking loop
for k = 1:maxIter
    
    %Calculates the hessian, if some eigenvalue of H is negative than
    %change the hessian for itself plus eye*|most negative eigenvalue + lambda|
    H = d2f(x(k,:));
    if sum(eig(H) < 0) > 0
        [L,D] = mcfac(H);
        H = L*D*L';
    end
    
    d = -inv(H)*df(x(k,:));  %Choose the direction as an steepest gradient descent
    alpha = 1;      %Step size, firs gess
    
    %Reduce the step until it satisfies the Armijo condition
    while f(x(k,:)' + alpha*d) > f(x(k,:)') + c1*alpha*df(x(k,:))'*d 
        alpha = ro*alpha;
    end
    
    %Find the next point
    x(k+1,:) = (x(k,:)' + alpha*d)';
end
out(:,n+1:2*n) = x
sol(:,2) = f(x(end,:)')
%% Backtraking Algorithm (Steepest Descent)

%Use this to find de derivative of a new function
% syms x;
% f = f(x);
% df = diff(f);

%Backtraking loop
for k = 1:maxIter
    d = -df(x(k,:));  %Choose the direction as an steepest gradient descent
    alpha = 1;      %Step size, firs gess
    
    %Reduce the step until it satisfies the Armijo condition
    while f(x(k,:)' + alpha*d) > f(x(k,:)') + c1*alpha*df(x(k,:))'*d 
        alpha = ro*alpha;
    end
    
    %Find the next point
    x(k+1,:) = (x(k,:)' + alpha*d)';
end
    
disp('Compare the steepest gradient descent vs newton method:')
out(:,2*n+1:3*n) = x
sol(:,3) = f(x(end,:)')
    
%Plot the error
error = zeros(maxIter,2);
for j = n:n:3*n
    for i = 1:maxIter
        error(i,j) = abs(f(out(i,j-n+1:j)') - mean(sol));
    end
end

plot(error(:,1), 'b'), hold on, plot(error(:,2), 'g'), plot(error(:,3), 'r')
legend('Newton Method With 2-norm','Newton Method With Cholesky','Steepest Gradient Descent')
title('Error of Steepest Gradient Descent vs Newton Method')
% axis([0 10 0 5])
    
    
    
    
    
    
    
    
    
    
    
    
    
    



