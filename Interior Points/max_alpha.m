function [alpha] = max_alpha(p, tau, s)
%Solve the problem:
%   max(alpha in (0,1): alpha*p >= -tau*s)
eps = 10^-4;        % tol

aux = -tau*s./p;
alpha = max(aux);
if alpha < 0
    disp('alpha < 0');
    alpha = eps;
elseif alpha > 1
    disp('alpha > 1');
    alpha = 1 - eps;
    
end

