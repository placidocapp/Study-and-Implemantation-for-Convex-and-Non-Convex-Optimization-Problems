function [x_opt,status] = solve_revised(A,b,c,x0,B,base,n,m)

while(1)
    u = B'\c(base);
    
    %Reduced costs
    s = c - A'*u;
    
    %Select next basic variable
    j = -1;
    for i = 1:n
        if s(i) < 0
            j = i;
        end
    end
    
    %Stop criteria
    if j == -1
        state = 'solved';
        break;
    end
    
    d = B\A(:,j);
    dp = d > 0;         %Positive values of d are indexed here
    
    %Verify if unbounded
    if sum(dp) == 0
        status = 'unbounded';
        break;
    end
    
    %Theta
    l = -1;
    min = 10^100;
    for i = 1:m
        if (dp(i) == 1) && (x(i)/d(i)) < min
            min = x(i)/d(i);
            l = i;
        end
    end
%     theta = min(x(dp)./d(dp));
    
    %Update x
    x = x - d'*theta;
    
    %Update B
    
    
    
end

end


