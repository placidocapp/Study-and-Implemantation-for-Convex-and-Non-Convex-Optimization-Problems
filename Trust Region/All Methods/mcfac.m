function [ L, D ] = mcfac( A )
%Modified Cholesky Factorization

%%  Parameters

eps = 10^-2;     %The factorization garantees that d >= eps
beta = 100;    %The factorization garantees that abs(m) <= beta

%% Inicialization
n = size(A,1);
if n ~= size(A,2) 
    disp('Error! The matrix A must be squared');
    return
end

L = eye(n);   %Triangular inferior
D = zeros(n);   %Diagonal matrix
c = zeros(n);   %Auxiliar matrix

%% Algorithm (Obtain LSL' = A)

for j = 1:n
    sum = 0;
    for s = 1:(j-1)
        sum = sum + D(s,s)*L(j,s)^2;
    end
    
    c(j,j) = A(j,j) - sum;
    D(j,j) = max([abs(c(j,j)), (max(max(abs(c)))/beta)^2, eps]);
    for i = (j+1):n
        sum = 0;
        for s = 1:(j-1)
            sum = sum + D(s,s)*L(j,s)*L(i,s);
        end
        c(i,j) = A(i,j) - sum;
        L(i,j) = c(i,j)/D(j,j);
    end
end

end

