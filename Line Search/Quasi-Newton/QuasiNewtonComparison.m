clear all;
close all;
clc

%%  Parameters

maxIter = 1000;  %Maximun number of iterations
maxSubIter = 30; %Max iterations of the bisection algorithm
n = 2;         %Dimension of the problem  
c1 = 10^-4;    %c1 and c2 are the constants of wolfe conditions    
c2 = 0.5;
eps = 10^-8;   %Stop criteria for the gradient
chooseAlgorithm = 1;    %If this flat equals to 0 uses de DFP method, 
                        %if 1 uses bfgs else uses regular newton step
chooseStepAlg = 2;      %If 0 choose the bisection else choose the zoom
                        %algorithm
                        
%% Inicializations
    
x = zeros(n,maxIter);   %Position of the search
% x(:,1) = randn(n,1);
x(:,1) = [-1.2 1];
p = zeros(n,maxIter);   %Direction 
H = zeros(n,n,maxIter); %Inverse of the aproxximate hessian
y = zeros(n,maxIter);   %yk = g(k+1) - g(k)
s = zeros(n,maxIter);   %Sk = x(k+1) - x(k) or sk = alphak*pk
kfinal = -1;
t = zeros(1,maxIter);   %Step lenght

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
for chooseAlgorithm = 0:2
    gkk = g(x(:,1));    %Gradient at initial point
    %Here i will use the hessian itself as initial point, but if you don't have
    %acess to it the initial gess should be approximately accurated, here is an
    %idea of how to do it
    % p(:,1) = -g(x(:,1));
    % alpha = 0.01/norm(p(:,1));
    % H(:,:,1) = pinv(alpha*p(:,1))'*( g(x(:,1)+alpha*p(:,1)) - g(x(:,1)) )';
    % H(:,:,1) = pinv(H(:,:,1));

    H(:,:,1) = pinv(B(x(:,1)));
    if sum(eig(H(:,:,1)) < 0) > 0
        [L,D] = mcfac(H(:,:,1));
        H(:,:,1) = L*D*L';
    end

%% Algorithm 

    kfinal = maxIter;
    for k = 1:maxIter
        %% Calculate the function value and it's derivatives
        fk = f(x(:,k));
        gk = gkk;

        %% Choose the direction
        p(:,k) = - H(:,:,k)*gk;

        %% Choose the step lenght by bisection method
        
        if chooseStepAlg == 0
            t(k) = bisec( f, g, x(:,k), p(:,k), c1, c2 );
        else
            t(k) = lineSearch( f, g, x(:,k), p(:,k), c1, c2 );
        end

        %% Stop criteria
        if norm(gk) < eps
            disp('Stopped at')
            kfinal = k
            break
        end

        %% Update

        x(:,k+1) = x(:,k) + t(k)*p(:,k);
        s(:,k) = t(k)*p(:,k);
        gkk = g(x(:,k+1));
        y(:,k) = gkk - gk;
        if chooseAlgorithm == 0
            %DFP Uptade
            H(:,:,k+1) = H(:,:,k) - (H(:,:,k)*y(:,k)*y(:,k)'*H(:,:,k))/...
                (y(:,k)'*H(:,:,k)*y(:,k)) + (s(:,k)*s(:,k)')/(y(:,k)'*s(:,k));
        elseif chooseAlgorithm == 1
            %BFGS Uptade
            ro = 1/(y(:,k)'*s(:,k));
            H(:,:,k+1) = ( eye(n) - ro*s(:,k)*y(:,k)' )...
                *H(:,:,k)*( eye(n) - ro*y(:,k)*s(:,k)' )...
                + ro*s(:,k)*s(:,k)';
        elseif chooseAlgorithm == 2
            %Newton method with modified hessian
            H(:,:,k+1) = pinv(B(x(:,k+1)));
            if sum(eig(H(:,:,k+1)) < 0) > 0
               [L,D] = mcfac(H(:,:,k+1));
               H(:,:,k+1) = L*D*L';
            end
        elseif chooseAlgorithm == 3
            %SR1 Update
            aux = (s(:,k)-H(:,:,k)*y(:,k));
            if abs( s'*aux ) >= 10^-8*norm(s)*norm(aux)   
                H(:,:,k+1) = H(:,:,k) + aux*aux'/( aux'*s(:,k) );
            else 
                H(:,:,k+1) = H(:,:,k);
            end
        end
    end

    x_final = x(:,kfinal)
    grad_norm = norm(g(x(:,kfinal)))
    f_final = f(x(:,kfinal))

    %% Plot some graphs

    x_error = zeros(kfinal,1);
    f_error = zeros(kfinal,1);
    g_evol = zeros(kfinal,1);
    for i = 1:kfinal
        x_error(i) = norm(sol) - norm(x(:,i));
        f_error(i) = f(sol') - f(x(:,i));
        g_evol(i) = norm(g(x(:,i)));
    end

    k = 1:kfinal;
   
    subplot(3,1,1), plot(k,x_error), hold on, title('Error in position x')
    legend('DFP','BFGS','Newton Method','SR1')
    subplot(3,1,2), plot(k,f_error), hold on, title('Error in function value')
    legend('DFP','BFGS','Newton Method','SR1')
    subplot(3,1,3), plot(k,g_evol), hold on, title('Gradient evolution')
    legend('DFP','BFGS','Newton Method','SR1')
end
    
    
    
    
    
    
    
    
    
    