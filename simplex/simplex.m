function [x_opt,f_opt,status] = simplex(c,A,b,Aeq,beq)
%simpex is na algorithm for solving linear problems
%Here A is already made of equality constrains
%       min c'*x
%       subject to
%           A*x = b

%% Variables

m = size(A,1) + size(Aeq,1);    %number of constrains
n = size(A,2);                  %number of variables
na = size(A,1);                 %number of slack variables
ny = m - size(A,1);             %number of extra variables

%Get together all restrictions
A = [A; Aeq];
b = [b; beq];

if m ~= length(b)
    disp('Number of elements in b and lines in A must be equal');
    return
end

if n ~= length(c)
    disp('Number of elements in c and columns in A must be equal');
    return
end
                                            
%% Fist Phase - Find an Feasible Solution

%If there are not enought slack variables, than use the extra variables 
%to find an initial solution
if ny > 0
    %Slack variables and aditional variables
    A = [A, eye(m)];
    c = [c; zeros(na,1)];
    n = n + ny + na;

    %Initial solution to auxiliary problem
    x0 = [zeros(n,1); b];

    %Vector with base variables
    base = (n+1):(n+m);
    
    %Auxiliary c for extra variables
    caux = c;
    c = [zeros(n-ny+m,1); ones(ny,1)];
    
    %Initial Tableau
    disp('Tableau inicial ...')
    tableau = [-c'*x0,    -c(base)'*A(:,1:(n-ny)), zeros(1,ny);
               x0(base),  A]
    
    %Solve the tableau
    [tableau, status] = solveTableau(tableau,m,n);
    
    disp('Tableau resolvido ...')
    tableau
    status
    
    if strcmp(status,'unbounded') == 1
        disp('It should not happen');
        x_opt = [];
        f_opt = -inf;
        return
    end
    
    %Function to identifý an base
    isbase = @(t,i) (( sum(t(:,i) == 0) == size(t(:,i),1)-1 ) && ...
            ( sum(t(:,i) == 1) == 1 ));
        
    %Garantee every extra variable is not in the base
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
                tableau = [tableau(1:(l-1),:); tableau((l+1):end,:)];
                m = m - 1;
                continue;
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
    
    disp('Tableau removendo variáveis extras ...')
    tableau
    
    %In the case we solved the tableau with extra variables correct than
    %just remove the extra variables
    tableau = tableau(:,1:(n+1-ny));
    
    %n returns to the original size
    n = n - ny;
    
    %The original c is now used
    c = caux;
    
    %Find x0, the initial solution
    x0 = zeros(n,1);
    base = zeros(m,1);
    cnt = 1;
    for i = 2:(n+1)
        if isbase(tableau,i)
            x0(i-1) = tableau(tableau(:,i) == 1,1);
            [aux, k] = max(tableau(:,i));   %insert in the base
            base(cnt) = k - 1;
            cnt = cnt + 1;
        end
    end
    
    %Now add the original cost to tableau
    tableau(1,1) = -c'*x0;
    for i = 2:(n+1)
        if ~isbase(tableau,i)
            tableau(1,i) = c(i-1) - c(base)'*tableau(2:end,i);
        end
    end
end
disp('Tableau entrando na fase 2 ...')
tableau

%% Phase 2

%If no extra variables are needed just start the tableau
if ny == 0
    tableau = [-c'*x0,    c';
               x0(base), A];
end

[tableau, status] = solveTableau(tableau,m,n);

disp('Tableau final ...')
tableau
%% Return

%Return x_opt and f_opt
x_opt = zeros(n-na,1);
for i = 2:(n+1-na)
    if isbase(tableau,i)
        x_opt(i-1) = tableau(tableau(:,i) == 1,1);
    end
end

f_opt = tableau(1,1);

end