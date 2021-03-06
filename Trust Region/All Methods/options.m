function [opt] = options(algorithm, hessAprox)
%Create a struct with the options for the algorithm
%algorithm     = 'cauchy', 'dogleg', 'subproblem'
%hessAprox     = 'none', 'sr1', 'dfp', 'bfgs'

%Default Structure
opt = struct(                              ...
    'Algorithm', 'dogleg',                 ...
    'Hessian_Aprox','none');

if exist('algorithm','var')
    opt.Algorithm = algorithm;
end

if exist('hessAprox','var')
    opt.Hessian_Aprox = hessAprox;

end

