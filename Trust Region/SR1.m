%This algorithm for SR1 can be find at pages 146 and 147 from the book
%Numerical Optimization. Nocedal, Jorge  and Wright, Stephen J. second
%edition

clear all;
close all;
clc

%% Parameters

maxIter = 100;       %Maximum number of iterations
n = 2;              %Dimension of the problem  
eps = 10^-4;        %Convergence tolerance
eta = 0.25;         %eta should belong to (0,1/4)
r = 10^0;           %r in (0,1)
delta0 = 1;         %Initial radius of trust region

%% Initializations

x = randn(n,maxIter);
Bk = zeros(n,n,maxIter);

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

% f = @(x) sin(x(1))+sin(x(2));
% g = @(x) [cos(x(1));
%           cos(x(2))];
% B = @(x) [  -sin(x(1)) 0;
%             0          -sin(x(2))];
% sol = [-1.570796327268957  -1.570796324699814];

% f = @(x) -200*exp(-0.2*sqrt(x(1)^2+x(2)^2));
% g = @(x) [(4*x(1).*exp(-(x(1).^2 + x(2).^2)^(1/2)/50))./(x(1).^2 + x(2).^2)^(1/2);
%           (4*x(2).*exp(-(x(1).^2 + x(2).^2)^(1/2)/50))./(x(1).^2 + x(2).^2)^(1/2)];
% B = @(x) [ (4*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(1/2) - (2*x(1)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2),                                                    - (2*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2);
%                                                     - (2*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2), (4*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(1/2) - (2*x(2)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(2)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2)] ;
% sol = [0 0];

%Rosenbrock
a = 1;
b = 100;
f = @(x) (a-x(1))^2 + b*(x(2)-x(1)^2)^2; 
g = @(x) [-2*(a - x(1))-4*b*x(1)*(x(2)-x(1)^2);
          2*b*(x(2)-x(1)^2)];
B = @(x) [2-4*b*x(2)+12*b*x(1)^2 -4*b*x(1);
          -4*b*x(1)                2*b];       
sol = [a a^2];

%% Initial gesses

gk = g(x(:,1));    %Gradient at initial point
%The book have a suggestion for inicialization in page 143, I prefer to
%calculate the hessian one time and use it as initial gess

Bk(:,:,1) = B(x(:,1));
if sum(eig(Bk(:,:,1)) < 0) > 0
    [L,D] = mcfac(Bk(:,:,1));
    Bk(:,:,1) = L*D*L';
end

%% Trust Region 

%Set the initial delta
delta = delta0;
kfinal = maxIter;
for k = 1:(maxIter-1)
    %% Dogleg
    
    %Here I use the dogleg method to make the method faster
    H = pinv(Bk(:,:,k));
    aux = gk'*Bk(:,:,k)*gk;
    pu = - gk'*gk/(aux)*gk;
    pb = - H*gk;
    
    %Find tau
    if norm(pb) <= delta
        %if pb is inside the trust region than we can go with it
        tau = 2;
    else
        %If tau is between 1 and 2 then we find the maximum feasible
        %tau
        %Closed formula solution
        aux = (pb - pu);
        a = aux'*aux;
        b = 2*pu'*aux;
        c = pu'*pu - delta^2;
        alpha = [(-b + sqrt(b^2-4*a*c))/(2*a);
                 (-b - sqrt(b^2-4*a*c))/(2*a)];
        if sum(imag(alpha)) > 0
            disp('numeric error')
            kfinal = k;
            break
        end
        tau = 1 + max(alpha);
    end
    
    %Choose next path direction based on tau 
    s = pu + (tau - 1)*( pb - pu );
    
    %% Algorithm
    
    y = g(x(:,k)+s) - gk;         %Calculathe y
    ared = f(x(:,k)) - f(x(:,k)+s);    %Actual reduction
    pred = -( gk'*s + 0.5*s'*Bk(:,:,k)*s );  %Predicted reduction
    ro = ared/pred;                    
    
    %Just update x if ro is positive
    if ro > eta
        x(:,k+1) = x(:,k) + s;
     
        %Update the gradient
        gk = g(x(:,k+1));
        if norm(gk) <= eps
            kfinal = k+1;
            break;
        end
    else
        x(:,k+1) = x(:,k);
    end
    
    %In case our confidance on the model is high we increase the area of
    %trust otherwise we decrease it
    if ro > 0.75
        if norm(s) > 0.8*delta
            delta = 2*delta;
        end
    elseif ro < 0.1 
        delta = 0.5*delta;
    end
   
    %Test the condition for update the inverse hessian
    aux = (y-Bk(:,:,k)*s);
    if abs( s'*aux ) >= r*norm(s)*norm(aux)   
        Bk(:,:,k+1) = Bk(:,:,k) + aux*aux'/( aux'*s + eps );   
    else
        Bk(:,:,k+1) = Bk(:,:,k);
    end
end

%Display the last values
x_algorithm = x(:,kfinal)
f_algorithm = f(x(:,kfinal))
x_sol = sol'
f_sol = f(sol')

%% Plot the solution

x_error = zeros(kfinal,1);
f_error = zeros(kfinal,1);
for i = 1:kfinal
    x_error(i) = norm(sol) - norm(x(:,i));
    f_error(i) = f(sol') - f(x(:,i));
end

k = 1:kfinal;
subplot(2,1,1), plot(k,x_error), title('Error in position x')
subplot(2,1,2), plot(k,f_error), title('Error in function value')        
        
        