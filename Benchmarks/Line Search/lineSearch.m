% MATLAB code by Plácido Campos based on Jorge Nocedal and Stephen J. Wright.
% Numerical Optimization, Second Edition.

function [x_opt, f_opt, time, gnorm, iter] = lineSearch(f,g,x0,B,opt,eps,maxIter,c1,c2,ro)
%The lineSearch function solve an optimization problem using the Line
%Search methods. To call lineSearch one should pass as arguments:
%           f: the function to optimize
%           g: the gradient of f
%           B: the Hessian of B, if not present an algorithm withought B
%              will be used
%           c1: first constant on wolfe conditions
%           c2: second constant on wolfe conditions
%           eps: expected norm of gradient
%           maxIter: max number of iterations
%           stepMethod:  'backtracking' or 
%                       'wolfe conditions'
%           ro: the constant for backtracking
%           modBMethod: Method for modifying the hessian. 
%                       'norm2' or 'modified cholesky'
tic
%% Define defatult variables

%If B is not defined it will not be used
if ~exist('B','var') || isempty(B)
    Bexist = 0;
else 
    Bexist = 1;
end

%The wolf conditions have some typical values
if ~exist('c1','var')
    c1 = 0.4;
end

if ~exist('c2','var')
    c2 = 10^-4;
end

%Default Options are 
if ~exist('opt','var')
    opt = options();
end

if ~exist('eps','var')
    eps = 10^-8;
end

if ~exist('maxIter','var')
    maxIter = 1000;
end

%If the used method is not backtracking this variable is left unused
if ~exist('ro','var')
    ro = 0.8;
end

%Newton need the B defined
if strcmp(opt.Algorithm,'newton') && Bexist == 0
    disp('Newton method needs the hessian...')
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
    end
end

if length(x0) ~= n
    disp('Gradient and x0 must have the same size')
end

%% Convert strings to flags

if strcmp(opt.Algorithm,'dfp')
    algorithm = 0;
elseif strcmp(opt.Algorithm,'bfgs')
    algorithm = 1;
elseif strcmp(opt.Algorithm,'newton')
    algorithm = 2;
elseif strcmp(opt.Algorithm,'sr1')
    algorithm = 3;
elseif strcmp(opt.Algorithm,'gradient')
    algorithm = 4;
else
    disp('Algorithm not recognized')
    return
end

if strcmp(opt.Step_Length_Method,'backtacking')
    stepMethod = 0;
elseif strcmp(opt.Step_Length_Method,'wolfeCond')
    stepMethod = 1;
elseif strcmp(opt.Step_Length_Method,'bisection')
     stepMethod = 2;
else   
    disp('Method for choose the step length not recognized')
    return
end

if strcmp(opt.Modified_Hessian_Method,'norm2')
    modMethod = 0;
elseif strcmp(opt.Modified_Hessian_Method,'cholesky')
    modMethod = 1; 
else   
    disp('Method for modify the hessian not recognized')
    return
end

%% Declarations
    
x = zeros(n,maxIter);   %Position of the search
x(:,1) = x0;
p = zeros(n,maxIter);   %Direction 
H = zeros(n,n,maxIter); %Inverse of the aproxximate hessian
y = zeros(n,maxIter);   %yk = g(k+1) - g(k)
s = zeros(n,maxIter);   %Sk = x(k+1) - x(k) or sk = alphak*pk
t = zeros(1,maxIter);   %Step lenght
kfinal = maxIter;       %The iteration that the algorithm stopped

%% Inicialization
    
gkk = g(x(:,1));    %Gradient at initial point

%If B is available calculate it the first time
if Bexist == 1
    H(:,:,1) = pinv(B(x(:,1)));
    if modMethod == 0
        if sum(eig(H(:,:,1)) < 0) > 0
            lambda = min(eig(H(:,:,1)));
            H(:,:,1) = H(:,:,1) + eye*(-lambda+0.1);
        end
    elseif modMethod == 1
        if sum(eig(H(:,:,1)) < 0) > 0
            [L,D] = mcfac(H(:,:,1));
            H(:,:,1) = L*D*L';
        end
    end
else
    %In case B is not available uses an gest
    p(:,1) = -gkk;
    alpha = 0.01/norm(p(:,1));
    H(:,:,1) = pinv(alpha*p(:,1))'*( g(x(:,1)+alpha*p(:,1)) - g(x(:,1)) )';
    H(:,:,1) = pinv(H(:,:,1));
end
    
%% Algorithm 

for k = 1:maxIter
    %% Calculate the function value and it's derivatives
    fk = f(x(:,k));
    gk = gkk;

    %% Choose the direction
    p(:,k) = - H(:,:,k)*gk;

    %% Choose the step lenght by bisection method

    if stepMethod == 0
        t(k) =  backtrack( f, g, x(:,k), p(:,k), c1, ro );
    elseif stepMethod == 1
        t(k) = wolfeCond( f, g, x(:,k), p(:,k), c1, c2 );
    elseif stepMethod == 2
        t(k) = bisec( f, g, x(:,k), p(:,k), c1, c2 );
    end

    %% Stop criteria
    if norm(gk) < eps
        disp('Stopped at')
        kfinal = k;
        break
    end

    %% Update

    x(:,k+1) = x(:,k) + t(k)*p(:,k);
    s(:,k) = t(k)*p(:,k);
    gkk = g(x(:,k+1));
    y(:,k) = gkk - gk;
    if algorithm == 0
        %DFP Uptade
        H(:,:,k+1) = H(:,:,k) - (H(:,:,k)*y(:,k)*y(:,k)'*H(:,:,k))/...
            (y(:,k)'*H(:,:,k)*y(:,k)) + (s(:,k)*s(:,k)')/(y(:,k)'*s(:,k));
    elseif algorithm == 1
        %BFGS Uptade
        ro = 1/(y(:,k)'*s(:,k));
        H(:,:,k+1) = ( eye(n) - ro*s(:,k)*y(:,k)' )...
            *H(:,:,k)*( eye(n) - ro*y(:,k)*s(:,k)' )...
            + ro*s(:,k)*s(:,k)';
    elseif algorithm == 2
        %Newton method with modified hessian
        H(:,:,k+1) = pinv(B(x(:,k+1)));
        
        if modMethod == 0
            if sum(eig(H(:,:,k+1)) < 0) > 0
                lambda = min(eig(H(:,:,k+1)));
                H(:,:,k+1) = H(:,:,k+1) + eye*(-lambda+0.1);
            end
        elseif modMethod == 1
            if sum(eig(H(:,:,k+1)) < 0) > 0
               [L,D] = mcfac(H(:,:,k+1));
               H(:,:,k+1) = L*D*L';
            end
        end
    elseif algorithm == 3
        %SR1 Update
        aux = (s(:,k)-H(:,:,k)*y(:,k));
        if abs( s'*aux ) >= 10^-8*norm(s)*norm(aux)   
            H(:,:,k+1) = H(:,:,k) + aux*aux'/( aux'*s(:,k) );
        else 
            H(:,:,k+1) = H(:,:,k);
        end
     elseif algorithm == 3
         %Gradient
         H(:,:,k+1) = eye(n);
    end
end

%% Return
x_opt = x(:,kfinal);
f_opt = f(x_opt);
disp(strcat(['Best objective found = ',num2str(f_opt)]))

disp(strcat(['Norm of gradient achieved = ',num2str(norm(g(x_opt)))]))

disp(strcat(['Number of iterations = ',num2str(kfinal)]))
time = toc;
gnorm = norm(g(x_opt));
iter = kfinal;
end

