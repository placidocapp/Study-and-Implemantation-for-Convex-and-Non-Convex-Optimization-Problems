function [tableau, status] = solveTableau(tableau,m,n)

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
        return;
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

end

