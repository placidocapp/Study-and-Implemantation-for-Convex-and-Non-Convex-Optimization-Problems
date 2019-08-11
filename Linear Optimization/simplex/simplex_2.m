function [x_opt,f_opt] = simplex_2(c,A,b,x0,base,H)
%simpex is na algorithm for solving linear problems
%Here A is already made of equality constrains
%       min c'*x
%       subject to
%           A*x = b

%% Variables

m = size(A,1);
n = size(A,2);

if m ~= length(b)
    disp('Number of elements in b and lines in A must be equal');
    return
end

if n ~= length(c)
    A
    c
    disp('Number of elements in c and columns in A must be equal');
    return
end

if ~exist('H','var')
    H = [];
    disp('falta completar o H')
    return
end

%% ALgorithm

x = x0;

for k = 1:1
    %Compute p
    p = (c(base)'*H)'
    
    %Reduced costs
    j = -1;
    for i = 1:(n)
        if(c(i)-p'*A(:,i) < 0)
            j = i
            break;
        end
    end
    
    %If all reduced costs are positive, we found the solution
    if j == -1 
        break
    end
    
    %Compute u
    u = H*A(:,j)
    
    %If there is not ui > 0 the solution ins -inf
    if max(u) <= 0 
        break
    end
    
    [theta, lb] = min(x(base)./u(u>0))
    l = base(lb)
    
    %Insert j in the base and remove l
    for i = 1:length(base)
        if base(i) == l
            base(i) = j;  
            break;
        end
    end
    base
    %Update x
    x(l) = 0;
    x(base) = x(base) - theta*u
    x(j) = theta
    
    A*x
    b
    f = c'*x
    
    %Update H
    H
    u
    H(lb,:) = H(lb,:)/u(lb);
    u(lb) = 1;
    for i = 1:length(u)
        if i ~= lb
            H(i,:) = H(i,:) - u(i);
            u(i) = 0;
        end
    end
    H
    u
    
end

x_opt = [];
f_opt = [];
end