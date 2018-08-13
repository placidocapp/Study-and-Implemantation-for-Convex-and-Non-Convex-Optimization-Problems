%This linear algorithms are based in the book Introduction to Linear
%Optimization from Bertsimas, Dimitris and Tsitsiklis, John N.
% Here I implement the Tableau Algorithm

clear all;
close all;
clc

%% Parameters


%% Inicialization


%% Problem definition 

%To begin with I will program only the constrains of the form Ax <= b so we
%have an feasible initial solution
c = -[10;12;12];
A = [1 2 2; 2 1 2; 2 2 1];
b = [20;20;20];
n = 3;          %Number of variables
nc = 3;         %Correction for n after add slack or surplus variables
m = size(b,1);  %Number of constrains

%Initial Tableau
T = [c'*zeros(n,1) [c' zeros(1,3)]; b  [A eye(m)]]

%% Full Tableau Algorithm

while(1)
    %First we will find the first negative reduced cost from the left for right
    %and pick the firs one we find to avoid cycling (Blend's rule)
    j = -1;
    for i = 2:(n+nc+1)
        if T(1,i) < 0
            j = i;
            break
        end
    end

    %Now j is the new basis
    if j == -1
        %Algorithm has found the best solution
        break;
    else
        %Algorithm continue, now we find the theta value, it will be the minor
        %value with positive value at the base
        minor_val = 10^10;  %The value of minor step size found
        minor_index = -1;   %The index of minor step size found
        for i = 2:(m+1)
            if T(i,j) > 0           %If the base is positive only
                aux = T(j,1)/T(j,i) %Step size
                if aux < minor_val  
                   minor_val = aux; %Update new best value
                   minor_index = i;
                end
            end
        end

        %If no theta is found so the cost is inf
        if minor_index == -1
            %Algorithm end and solution is inf or -inf
            disp('The sol in inf of -inf')
            break;
        else
            theta = minor_val;
            i = minor_index;
        end
        i
        j
        %Now we know that (i,j) is going to enter the new basis, we just
        %need to turn the new base
        T = changebase( T, i, j );
        T
                  
    end
    pause
end

%% Change the basis

%In some cases the final tableau will have some of the slack variables at
%the basic solution, in this case we just change tha base

for k = nc+1:n+nc+1
    if isbase(T(:,k))
        %Choose the best new base
        min_val = 10^10;
        j = -1;
        for p = 2:n+1
            if T(1,p) < min_val && isbase(T(:,p)) == 0
                min_val = T(1,p);
                j = p;
            end
        end
        
        %Since we know column k is base the bigger term is 1
        [max_val, i] = max(T(2:end,k));
        
        %Now we know the new column and line j and i
        T = changebase( T, i+1, j );
        T
    end
end



    
    
    
    
    
    
    
    



