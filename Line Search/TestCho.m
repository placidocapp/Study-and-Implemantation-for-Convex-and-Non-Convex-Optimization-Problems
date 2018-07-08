clear all
close all
clc

eps = 10^-10;
for i = 1:1000
    A = 100*randn(5);
    A = A'*A;
    A = A+A';
    [L,D] = mcfac(A);
    C = L*D*L';
    eig(A)
    eig(C)
    for k = size(A,1)
        for j = size(A,2)
            if (A(k,j) > C(k,j)+eps) || ((A(k,j) < C(k,j)-eps)) 
                disp('erro');
                return
            end
        end
    end
    v = eig(C);
    for k = 1:size(v,1)
        if (v(k) < 10^-2)
            disp('Def negativa, erro !')
            return
            A
        end
    end
    i
end