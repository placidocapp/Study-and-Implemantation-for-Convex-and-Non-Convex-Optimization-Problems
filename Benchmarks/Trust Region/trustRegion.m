function [x_opt, f_opt, time, gnorm, iter] = trustRegion(f,g,x0,B,opt,eps,maxIter,eta,delta0...
                                     ,r, maxIterSub)
%The trustRegion function solve an optimization problem using the Trust
%Regiom methods. To call trustRegion one should pass as arguments:
%           f: the function to optimize
%           g: the gradient of f
%           x0: start point
%           B: the Hessian of B, if not present an algorithm withought B
%              will be used
%           opt: select the methods options
%           eps: expected norm of gradient
%           maxIter: max number of iterations
%           eta: preferably between (0,1/4), if ro < eta than the region is
%           not trustable
%           delta0: the initial gess for the trust area
%           r: preferably between (0,1)
%           maxIterSub: Max iterations to the subproblem method
tic
%% Define defatult variables

%If B is not defined it will not be used
if ~exist('B','var') || isempty(B)
    Bexist = 0;
else 
    Bexist = 1;
end

%The wolf conditions have some typical values
if ~exist('eta','var')
    eta = 0.25;
end

if ~exist('delta0','var')
    delta0 = 1;
end

%Default Options are 
if ~exist('opt','var')
    opt = options();
end

if ~exist('eps','var')
    eps = 10^-20;
end

if ~exist('maxIter','var')
    maxIter = 1000;
end

if ~exist('maxIterSub','var')
    maxIterSub = 10;
end

%If the used method is not backtracking this variable is left unused
if ~exist('r','var')
    r = 1;
end

%Newton need the B defined
if strcmp(opt.Hessian_Aprox,'none') && Bexist == 0
    disp('Hessian is needed...')
    x_opt = [];
    f_opt = [];
    return
end

%% Define the size of the problem
n = length(x0);

%Test size
if Bexist == 1
    if size(B(x0),1) ~= size(B(x0),2) && size(B(x0),1) ~= n
        disp('Hessian sizes must be equal to columns in gradient')
        return
    end
end

if length(x0) ~= n
    disp('Gradient and x0 must have the same size')
    return
end
%% Convert strings to flags

if strcmp(opt.Algorithm,'cauchy')
    algorithm = 0;
elseif strcmp(opt.Algorithm,'dogleg')
    algorithm = 1;
elseif strcmp(opt.Algorithm,'subproblem')
    algorithm = 2;
else
    disp('Algorithm not recognized')
    return
end

if strcmp(opt.Hessian_Aprox,'none')
    hessAprox = 0;
elseif strcmp(opt.Hessian_Aprox,'sr1')
    hessAprox = 1;
elseif strcmp(opt.Hessian_Aprox,'dfp')
    hessAprox = 2;
elseif strcmp(opt.Hessian_Aprox,'bfgs')
    hessAprox = 3;
else   
    disp('Method for approximate the hessian not recognized')
    return
end
%% Declarations
    
x = randn(n,maxIter);
x(:,1) = x0;
Bk = zeros(n,n,maxIter);
maxIterSub = 30;
lambda = zeros(maxIterSub,1);  %Trust Region subproblem
lambda(1) = 1;

%% Inicialization
    
gk = g(x(:,1));    %Gradient at initial point

%If B is available calculate it the first time
if Bexist == 1
    Bk(:,:,1) = B(x(:,1));
    
    %Modified Cholesky
    if sum(eig(Bk(:,:,1)) < 0) > 0
        [L,D] = mcfac(Bk(:,:,1));
        Bk(:,:,1) = L*D*L';
    end
else
    %In case B is not available uses an gest
    p = -gk;
    alpha = 0.01/norm(p);
    Bk(:,:,1) = pinv(alpha*p)'*( g(x(:,1)+alpha*p) - g(x(:,1)) )';
end
%% Trust Region 

%Set the initial delta
delta = delta0;
kfinal = maxIter;
for k = 1:(maxIter-1)
    if algorithm == 0
        %Cauchy algorithm
        aux = gk'*Bk(:,:,k)*gk;
        if aux <= 0
            tau = 1;
        else 
            tau = min( 1, norm(gk)^3/(delta*aux) );
        end
        s = - tau*delta*gk/norm(gk);
    
    elseif algorithm == 1
        % Dogleg
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
            alpha = [(-b + sqrt(b^2-4*a*c))/(2*a),...
                     (-b - sqrt(b^2-4*a*c))/(2*a)];
            if sum(imag(alpha) ~= 0) > 0
                disp('numeric error')
                kfinal = maxIter;
                break
            end
            tau = 1 + max(alpha);
        end

        %Choose next path direction based on tau 
        s = pu + (tau - 1)*( pb - pu );
    elseif algorithm == 2
        %Trust Region Subproblem
        [L,D] = mcfac(Bk(:,:,k));
        Bk(:,:,k) = L*D*L';
        
        for i = 1:maxIterSub
            
            %If B + lambda*I <= 0 then correct lambda
            aux2 = eig(Bk(:,:,k)+lambda(i)*eye(n));
            if sum(aux2 < 0) > 0
                lambda(i) = lambda(i)-(1+10^-3)*min(aux2);
            end
            
            aux2 = eig(Bk(:,:,k)+lambda(i)*eye(n));
            if sum(aux2 < 0) > 0
                x_opt = x(:,k);
                f_opt = f(x(:,k));
                time = toc;
                gnorm = norm(g(x_opt));
                iter = maxIter;
                return
            end

            %Call cholesky decomposition
            R = chol(Bk(:,:,k)+lambda(i)*eye(n));
            s = -pinv(R'*R)*gk;
            q = pinv(R')*s;
            lambda(i+1) = lambda(i) + ( norm(s)/norm(q) )^2*...
            ( ( norm(s)-delta )/delta ); 
            
            %Stop criteria
            if abs(lambda(i+1) - lambda(i)) < eps
                break
            end
            
            %Numeric error
            if abs(lambda(i+1)) > 10^300
                disp('numeric error')
                kfinal = maxIter;
                break 
            end
        end
        
        if kfinal ~= maxIter
            kfinal = maxIter;
            break;
        end
       
        lambda(1) = 1;
    end
    %% Hessian Approximation
    
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
    
    %Stop criteria
    if norm(gk) < eps
        kfinal = k;
        break
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
   
    if hessAprox == 0
        Bk(:,:,k+1) = B(x(:,k));
        [L,D] = mcfac(Bk(:,:,k));
        Bk(:,:,k) = L*D*L';
    elseif hessAprox == 1
        %SR1
        %Test the condition for update the inverse hessian
        aux = (y-Bk(:,:,k)*s);
        if abs( s'*aux ) >= r*norm(s)*norm(aux)   
            Bk(:,:,k+1) = Bk(:,:,k) + aux*aux'/( aux'*s + eps );   
        else
            Bk(:,:,k+1) = Bk(:,:,k);
        end
    elseif hessAprox == 2
        %DFP
        if norm(s) > 10^-100
            rok = 1/( y'*s );
            Bk(:,:,k+1) = (eye(n)-rok*y*s')*Bk(:,:,k)*(eye(n)-rok*s*y')...
                            + rok*(y*y');
        else
            Bk(:,:,k+1) = Bk(:,:,k);
        end
    elseif hessAprox == 3
        %BFGS
        %if the trust region size is changing
        if norm(s) > 10^-100
            Bk(:,:,k+1) = Bk(:,:,k) - ...
                (Bk(:,:,k)*(s*s')*Bk(:,:,k))/(s'*Bk(:,:,k)*s) + ...
                (y*y')/(y'*s);
        else 
             Bk(:,:,k+1) = Bk(:,:,k);
        end
    end
    
    %Nan Case
    if (sum(sum(Bk(:,:,k+1))) == inf) || (isnan(sum(sum(Bk(:,:,k+1)))) == 1)
        x_opt = x(:,k);
        f_opt = f(x(:,k));
        time = toc;
        gnorm = norm(g(x_opt));
        iter = maxIter;
        return;
    end
        
end

%% Return
x_opt = x(:,kfinal);
f_opt = f(x_opt);
disp(strcat(['Best objective found = ',num2str(f_opt)]))

disp(strcat(['Norm of gradient achieved = ',num2str(norm(g(x_opt)))]))

disp(strcat(['Number of iterations = ',num2str(kfinal)]))

if norm(g(x_opt)) < eps
    disp('Status: Solved')
elseif norm(g(x_opt))/100 < eps
    disp('Status: Badly Solved')
else 
    disp('Status: Not Solved, bad convergence')
end
    
time = toc;
gnorm = norm(g(x_opt));
iter = kfinal;


end

