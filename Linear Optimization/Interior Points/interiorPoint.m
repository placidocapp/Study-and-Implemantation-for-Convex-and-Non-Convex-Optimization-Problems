function [x_opt,f_opt,status] = interiorPoint(c,A,b,Aeq,beq)
%This interior Point Method solves a problem of the form
%       min c'*x
%       subject to
%           A*x <= b
%           Aeqx = beq

%% Variables

m = size(A,1) + size(Aeq,1);    %number of constrains
n = size(A,2) + size(Aeq,2);    %number of variables
na = size(A,1);                 %number of slack variables

if m ~= length(b) + length(beq)  && ~isempty(A)
    disp('Number of elements in b and lines in A must be equal');
    return
end

if n ~= length(c) + na && ~isempty(A)
    disp('Number of elements in c and columns in A must be equal');
    return
end

%Get together all restrictions
A = [A; Aeq];
b = [b; beq];

%Add slack variables
c = [c; zeros(na,1)];
A = [A eye(m,na)];

%Precision
eps = 10^-4;

%Max Iterations
maxIter = 1000;

%% Initial Solution (Heuristic)

%Initial values, satisfy the restrictions but not the 
aux = pinv(A*A');
x = A'*aux*b;
l = aux*A*c;
s = c - A'*l;

%Now values of x and s will be positive
x = x + max( -1.5*min(x), 0);
s = s + max( -1.5*min(s), 0);

%Garantee x0 and s0 are not too close to zero and not too dissimilar
x = x + ( 0.5*(x'*s)/sum(s) );
s = s + ( 0.5*(x'*s)/sum(x) );

%% Algorithm

%First mi
mi = x'*s;

%First eta
eta = 0.9;

%First Jacobian
J = [zeros(n,n), A', eye(n,n);
     A, zeros(m,m), zeros(m,n);
     diag(s), zeros(n,m), diag(x)];

final_k = maxIter;
for k = 1:maxIter
    
    %Jacobian Update
    J((end-n+1):end,:) = [diag(s), zeros(n,m), diag(x)];
     
    %Get affine delta
    F = [ c-s-A'*l; b-A*x; -diag(x)*diag(s)*ones(n,1)];
%     da = J\F;
%     da = pinv(J)*F;
    da = linsolve(J,F);

    %Calculate the step length
    aux = da(1:n,1) < 0;
    ax = max(min(1, min(-x(aux)./da(aux,1))),0.4);
    aux1 = da((end-n+1):end,1) < 0;
    aux = [zeros(n+m,1) == 1; aux1];
    as = max(min(1, min(-x(aux1)./da(aux,1))),0.4);
    if isempty(as) || isempty(ax)
        status = 'not solved';
        x_opt = x;
        f_opt = -inf;
        return
    end
    
    %Value of mi that will be obtained
    miaff = (x + ax*da(1:n,1))'*(s + as*da((end-n+1):end,1))/n;
    
    %Centering parameter
    sigma = (miaff/mi)^3;
    mi = x'*s;
    
    %Stop Criteria
    if mi < eps
        status = 'solved';
        final_k = k;
        break;
    end
    
    %Find the direction
    F((end-n+1):end, :) = F((end-n+1):end, :) + sigma*mi -...
                   diag(da(1:n,:))*diag(da((end-n+1):end,:))*ones(n,1);
%     d = J\F;
    d = linsolve(J,F);
    
    %Calculate the new step lengths
    aux = d(1:n,1) < 0;
    ax = max(min(1, min(-x(aux)./d(aux,1))),0.4);
    aux1 = d((end-n+1):end,1) < 0;
    aux = [zeros(n+m,1) == 1; aux1];
    as = max(min(1, min(-x(aux1)./d(aux,1))),0.4);
    if isempty(as) || isempty(ax)
        status = 'not solved';
        x_opt = x;
        f_opt = -inf;
        return
    end
    
    %Finaly the step lengths
    eta = min(eta*1.01, 1);
    ax = min(1, eta*ax);
    as = min(1, eta*as);

    %Update the position
    x = x + ax*d(1:n,:);
    s = s + as*d((end-n+1):end,:);
    l = l + as*d((n+1):(n+m),:);
end

%Return the solution
x_opt = x;
f_opt = c'*x;

if final_k == maxIter
    status = 'bad convergence';
    return
end

if sum(A*x ~= b) > 0
    status = 'infeasible';
end

end

