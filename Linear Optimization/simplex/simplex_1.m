function [x_opt,f_opt] = simplex_1(c,A,b,x0,base,B)
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

if ~exist('B','var')
    B = [];
    disp('falta completar o B')
    return
end

%% ALgorithm

x = x0;

for k = 1:2
    %Compute p
    p = (c(base)\B)'
    
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
    u = B\A(:,j)
    
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
    x(base) = x(base) - theta*u
    x(j) = theta
    A
    A*x
    b
    
    %Update H
    B(:,lb) = A(:,j)
    
    f = c'*x
    
end

x_opt = [];
f_opt = [];
end