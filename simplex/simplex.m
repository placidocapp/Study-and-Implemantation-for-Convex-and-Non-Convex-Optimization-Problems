function [x_opt,f_opt,status] = simplex(c,A,b,x0,base,ny)
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

%Number of extra variables
if ~exist('ny','var')
    ny = 0; 
end

if ny > 0
    tableau = [-c'*x0,    -c(base)'*A(:,1:(n-ny)), zeros(1,ny);
           x0(base), A]
else
    tableau = [-c'*x0,    c';
               x0(base), A];
end

%% ALgorithm

while(1)
    %Find an negative reduced cost
    j = -1;
    for i = 2:(n+1)
        if tableau(1,i) < 0
            j = i;
            break;
        end
    end
    j;
    %If the reduced costs are nonegative we found the solution
    if j == -1
        status = 'solved';
        break;
    end
    
    %Find theta
    min = 10^100;
    for i = 2:(m+1)
        if tableau(i,j) > 0 && tableau(i, 1)/tableau(i,j) < min
            min = tableau(i, 1)/tableau(i,j);
            l = i;
        end
    end
    l;
    %It there is not an positive component of u, the solution is 
    %unbounded
    if min == 10^100;
        status = 'unbounded';
        x_opt = [];
        f_opt = -inf;
        return
    end
    
    %Update the tableau
    tableau(l,:) = tableau(l,:)/tableau(l,j);
    for i = 1:(m+1)
        if i ~= l
            tableau(i,:) = tableau(i,:) - ...
                tableau(l,:)*tableau(i,j)/tableau(l,j);
        end
    end
    tableau;
end

%Function to identif� an base
isbase = @(t,i) (( sum(t(:,i) == 0) == size(t(:,i),1)-1 ) && ...
            ( sum(t(:,i) == 1) == 1 ));

%if this is phase 1, garantee every extra variable is not in the base
if ny ~= 0
    tableau
    for i = (n+2-ny):(n+1)
        %if it's base than try changing the base
        if isbase(tableau,i)
            i;
            %Find line for switch 
            for k = 2:(m+1)
                if tableau(k,i) == 1
                    l = k;
                    break;
                end
            end
            l;
            
            %try to find an reduced cost = 0
            j = -1;
            for k = 2:(n+1-ny)
                if tableau(1,k) == 0 && ~isbase(tableau,k)
                    j = k;
                    break;
                end
            end
            j;
            %If the solution needs the extra variables
            if j == -1
                status = 'infeasible';
                x_opt = [];
                f_opt = [];
                return
            end
            
            %If tableau(l,j) = 0 than we have lines LD. Remove line l
            if tableau(l,j) == 0
                tableau = [tableau(1:(l-1),:); tableau((l+1):end,:)]
                break;
            end
            
            %if it's possible, change the base
            tableau(l,:) = tableau(l,:)/tableau(l,j);
            for k = 1:(m+1)
                if k ~= l
                    tableau(k,:) = tableau(k,:) - ...
                        tableau(l,:)*tableau(k,j)/tableau(l,j);
                end
            end
            tableau;
        end
    end
end

%Return x_opt and f_opt
xopt = zeros(n,1);
for i = 2:(n+1)
    if isbase(tableau,i)
        x_opt(i-1) = tableau(tableau(:,i) == 1,1);
    end
end

f_opt = tableau(1,1);

end