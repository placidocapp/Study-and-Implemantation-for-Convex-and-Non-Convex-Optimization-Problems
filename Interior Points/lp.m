function [x_opt,f_opt,status] = lp(c,A,b,Aeq,beq)
%lp is a function for solving linear problems
%The entrys are similar to linprog (matlab function)
%           f: the function to me minimized
%           A and b: The values to form an inequality system Ax <= b
%           Aeq and beq: The values to form an equality system Aeqx <= beq
tic
%% Define variables
if ~exist('Aeq','var')
    Aeq = [];
end
   
if ~exist('beq','var')
    beq = [];
end

if ( ~isempty(A) && isempty(b) ) || ( isempty(A) && ~isempty(b) )
    disp('If A exist, b must exist too...')
    return
end

if ( ~isempty(Aeq) && isempty(beq) ) || ( isempty(Aeq) && ~isempty(beq) )
    disp('If Aeq exist, beq must exist too...')
    return
end

if size(Aeq,2) ~=  + size(A,2) && ~isempty(Aeq) && ~isempty(A)
    disp('Number of variables in Aeq must be equal to A...')
    return
end

%% Call Simplex
[x_opt,f_opt,status] = interiorPoint(c,A,b,Aeq,beq);
toc
end

