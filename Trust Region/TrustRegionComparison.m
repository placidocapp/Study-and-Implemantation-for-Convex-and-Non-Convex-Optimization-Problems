% Algorithm from Numerical Optimization. Jorge Nocedal and Stephen J. Wright.
% Second edition. This algorithm can be found at pages 69, and the way to
% choose pk can be found at pages 71-77
% MATLAB code by Plácido Campos

clear all;
close all;
clc

%%  Parameters

maxIter = 50;  %Maximun number of iterations
n = 2;          %Dimension of the problem   
eta = 0.20;     %Parameter that decides if the algorithm steps or not,
                %the book recomends it between (0,0.25)
how_choose_step = 1; %If 0 then choose the step with cauchy, else use 
                     %dogleg
maxIterSub = 10;    %Max Iterations of the subproblem
eps = 10^-8;        %Stop criteria to Trust region subptoblem

%% Inicialization

p = zeros(maxIter,n);       %Direction of the step
ro = zeros(maxIter,1);      %rate of real decrease divided by estimatede decrease
x = zeros(maxIter,n);       %Actual point
x(1,:) = randn(1,n);
delta = zeros(maxIter,1);   %Radius of trust area
delta(1) = 1;
deltaMax = 4;
tau = zeros(maxIter,1);
B = zeros(n,n,maxIter);     %Hessian or modified hessian
k = 1;
pu = zeros(maxIter,n);      %Dogleg mehtod
pb = zeros(maxIter,n);      %Dogleg mehtod
syms tau_aux;               %Dogleg method
lambda = zeros(maxIterSub,1);  %Trust Region subproblem
lambda(1) = 1;
q = zeros(n,1);             %Trust Region subproblem

%% Function

%Function values
A = randn(n);
A = A'*A;
A = A'+ A;
b = randn(n,1);

%Calculate the function value and its derivatives
% f = @(x) 0.5*x'*A*x + b'*x;
% g = @(x) A*x + b;
% H = @(x) A;
% sol = -(pinv(A)*b)';

%Trigonometric function
% f = @(x) sin(x(1))+sin(x(2));
% g = @(x) [cos(x(1));
%           cos(x(2))];
% H = @(x) [  -sin(x(1)) 0;
%             0          -sin(x(2))];
% sol = [-1.570796327268957  -1.570796324699814];

% f = @(x) -200*exp(-0.2*sqrt(x(1)^2+x(2)^2));
% g = @(x) [(4*x(1).*exp(-(x(1).^2 + x(2).^2)^(1/2)/50))./(x(1).^2 + x(2).^2)^(1/2);
%           (4*x(2).*exp(-(x(1).^2 + x(2).^2)^(1/2)/50))./(x(1).^2 + x(2).^2)^(1/2)];
% H = @(x) [ (4*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(1/2) - (2*x(1)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2),                                                    - (2*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2);
%                                                     - (2*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(1)*x(2)*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2), (4*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(1/2) - (2*x(2)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(25*(x(1)^2 + x(2)^2)) - (4*x(2)^2*exp(-(x(1)^2 + x(2)^2)^(1/2)/50))/(x(1)^2 + x(2)^2)^(3/2)] ;
% sol = [0 0];


%Rosenbrock
a = 1;
b = 100;
f = @(x) (a-x(1))^2 + b*(x(2)-x(1)^2)^2; 
g = @(x) [-2*(a - x(1))-4*b*x(1)*(x(2)-x(1)^2);
          2*b*(x(2)-x(1)^2)];
H = @(x) [2-4*b*x(2)+12*b*x(1)^2 -4*b*x(1);
          -4*b*x(1)                2*b];       
sol = [a a^2];

m = @(x,p) f(x) + g(x)'*p + 0.5*p'*B(:,:,k)*p;

%% Algorithm
for how_choose_step = 0:2
    for k = 1:maxIter
        %% Calculate the values of gradient and hessian
        B(:,:,k) = H(x(k,:)');
        gk = g(x(k,:)');
        fk = f(x(k,:)');

        %Check if B is semidefinite positive, if not then correct it with
        %modified cholesky factorization
        if sum(eig(B(:,:,k)) < 0) > 0
            [L,D] = mcfac(B(:,:,k));
            B(:,:,k) = L*D*L';
        end

        %% Solve argmin mk(p) to find pk
        aux = gk'*B(:,:,k)*gk;
        %Cauchy algorithm
        if how_choose_step == 0
            if aux <= 0
                tau(k) = 1;
            else 
                tau(k) = min( 1, norm(gk)^3/(delta(k)*aux) );
            end
            p(k,:) = - tau(k)*delta(k)*gk/norm(gk);
        elseif how_choose_step == 1
            %Dogleg Method from Numerical optimization page 74
            pu(k,:) = - gk'*gk/(aux)*gk;
            pb(k,:) = - pinv(B(:,:,k))*gk;

            %Find tau
            if norm(pb(k,:)) <= delta(k)
                %if pb is inside the trust region than we can go with it
                tau(k) = 2;
            else
                %If tau is between 1 and 2 then we find the maximum feasible
                %tau
                %Closed formula solution
                aux = (pb(k,:) - pu(k,:))';
                a = aux'*aux;
                b = 2*pu(k,:)*aux;
                c = pu(k,:)*pu(k,:)' - delta(k)^2;
                alpha = [(-b + sqrt(b^2-4*a*c))/(2*a);
                         (-b - sqrt(b^2-4*a*c))/(2*a)];
                if sum(imag(alpha)) > 0
                    disp('numeric error')
                    kfinal = k;
                    break
                end
                tau(k) = 1 + max(alpha);
            end

            %Choose next path direction based on tau
            p(k,:) = pu(k,:) + (tau(k) - 1)*( pb(k,:) - pu(k,:) );
        else
            %Trust Region Subproblem
            for i = 1:maxIterSub
               %If B + lambda*I <= 0 than correct lambda
                aux_ = eig(B(:,:,k)+lambda(i)*eye(n));
                if sum(aux_ < 0) > 0
                    lambda(i) = lambda(i)-min(aux_)-eps*min(aux_);
                end

                R = chol(B(:,:,k)+lambda(i)*eye(n));
                p(k,:) = -pinv(R'*R)*gk;
                q = pinv(R')*p(k,:)';
                lambda(i+1) = lambda(i) + ( norm(p(k,:))/norm(q) )^2*...
                ( ( norm(p(k,:))-delta(k) )/delta(k) ); 

                %Stop creteria
                if abs(lambda(i+1) - lambda(i)) < eps
                    break
                end
            end
            %Gess the initial lambda equals to the last one
            lambda(1) = lambda(i);
        end
        %% Decide to step or to rise/decrease the area delta

        %Evaluate ro(k)
        ro(k) = ( fk - f(x(k,:)'+p(k,:)') )/...
            ( m(x(k,:)',zeros(n,1)) - m(x(k,:)',p(k,:)') );
        if ro(k) < 0.25
            delta(k+1) = 0.25*delta(k);
        else
            if (ro(k) > 3/4 || norm(p(k)) == delta(k))
                delta(k+1) = min(2*delta(k), deltaMax);
            else
                delta(k+1) = delta(k);
            end
        end

        %Decide to take the step or not
        if ro(k) > eta
            x(k+1,:) = x(k,:) + p(k,:);
        else
            x(k+1,:) = x(k,:);
        end
    end

    x(end,:)
    f(x(end,:)')

    %% Plot some graphs
    x_error = zeros(maxIter,1);
    f_error = zeros(maxIter,1);
    for i = 1:maxIter
        x_error(i) = norm(sol) - norm(x(i,:));
        f_error(i) = f(sol') - f(x(i,:)');
    end

    k = 1:maxIter;
    subplot(2,1,1), plot(k,x_error), title('Error in position x'),hold on
    legend('Cauchy Point','Dogleg','Subproblem')
    subplot(2,1,2), plot(k,f_error), title('Error in function value'),hold on
    legend('Cauchy Point','Dogleg','Subproblem')
end
        
        
        
        
        
        
        
        
        
        