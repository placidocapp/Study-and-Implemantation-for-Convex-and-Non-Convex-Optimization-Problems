clear all;
close all;
clc

%% Problem

%   min f(x)
%   subject to
%       h(x)     =  0
%       g(x) + s =  0 -> g(x) <= 0
%              s >= 0 

% Equality constrais --------------------
% Ah = [2 1;
%      1 2];
% bh = [1; 1];
% h = @(x) x'*Ah*x + bh'*x;
% hx = @(x) (2*Ah*x + bh)';
% hxx = @(x) 2*Ah;

h = @(x) 0;
hx = @(x) [0;0]';
hxx = @(x) [0 0;0 0];


% Inequality constrains -----------------
% Ag = [2 1;
%      1 2];
% bg = [1; 1];
% g = @(x) x'*Ag*x + bg'*x;
% gx = @(x) (2*Ag*x + bg)';
% gxx = @(x) 2*Ag;

% g = @(x) [x(1); x(2)];       % x1 <= 0 ; x2 <= 0
% gx = @(x) [1 0; 0 1];
% gxx = @(x) cat(3,[0 0; 0 0],[0 0; 0 0]);

g = @(x) 0;       % x1 <= 0 ; x2 <= 0
gx = @(x) [0;0]';
gxx = @(x) cat(3,[0 0; 0 0],[0 0; 0 0]);

% Objective function --------------------
Af = [2 1;
     1 2];
bf = [1; 1];
f = @(x) x'*Af*x + bf'*x;
fx = @(x) 2*Af*x + bf;
fxx = @(x) 2*Af;

% Initial point --------------------
x0 = [0.5; -0.5]

%% Call solver

n = length(x0);
mh = length(h(x0));
mg = length(g(x0)); 

%% Parameters

mi = 0.1;
eta = 0.5;
sigma = 0.5;
eps_mi = 0.01;
eps_tol = 0.01;
tau = 0.995;
p = 1;          % Norm to use, select 1 or 2
ro = 0.8;       % Backtrack parameter
po = 0.9;       % Value for barrier parameter
c1 = 10^-2;     % Parameter for armijo condition  

%% Initial Solution

s0 = rand(mg,1);
S = diag(s0);
e = ones(mg,1);

% Calculate z0
z0 = mi*pinv(S)*e;
y0 = hx(x0)'\( gx(x0)'*z0 - fx(x0) );

% Initializations
delta = 0;
x = x0;
s = s0;
y = y0';
z = z0;

while E(x,s,y,z,0,p,fx,hx,gx,h,g) > eps_tol
    while E(x,s,y,z,mi,p,fx,hx,gx,h,g) > eps_mi
        % C�lculo de p --------------------------
        L = fxx(x) - mul(y,hxx(x)) - mul(z,gxx(x));
        
        F = [L + delta*eye(n) zeros(n,mg) hx(x)' gx(x)';
             zeros(mg,n) pinv(diag(s))*diag(z) zeros(mg,mh)   -eye(mg);
             hx(x)   zeros(mh,2*mg+mh);      
             gx(x)  -eye(mg)  zeros(mg,mh+mg)];
        
        % Conditioning F matrix -----------------
        F_eig = eig(F);
        % Verify the inertia of F, we must have 0 0's eigenvalues and n +
        % mg positive eigenvalues
        if sum(F_eig == 0) ~= 0 || sum(F_eig >= 0) ~= n + mg
            % Add si for conditioning
            if sum(F_eig == 0) > 0
                si = 10^-8*eta*mi^beta
                F(n+mg+1:n+mg+mh,n+mg+1:n+mg+mh) = eye(mh)*si
            end
            
            % Initialize delta
            if delta == 0
                delta = 10^-4;
            else
                delta = delta/2;
            end
            
            %Increase delta until we find the required inertia for F
            while 1
                F_eig = eig(F)
                if sum(F_eig == 0) == 0 && sum(F_eig >= 0) == n + mg
                    break
                else
                    delta = delta*10;
                end
                
                % Update F
                F(1:n,1:n) = L + delta*eye(n);
            end
        end
        
        % Calculating directions p --------------
        D = [fx(x) - hx(x)'*y - gx(x)'*z;
             z - mi*pinv(diag(s))*ones(mg,1);
             h(x);
             g(x) - s];
        P = F\D;
        px = P(1:n);
        ps = P(n+1:n+mg);
        py = P(n+mg+1:n+mg+mh);
        pz = P(n+mg+mh+1:n+2*mg+mh);
        
        % Calculating step size -----------------
        alphasmax = max_alpha(ps, tau, s);
        alphazmax = max_alpha(pz, tau, z);
        
        alphas = backtrack(x,s,px,ps,ro,po,p,alphasmax,f,fx,h,hx,g,gx,mi,c1);
        alphaz = alphazmax;
        
        % Update --------------------------------
        x = x + alphas*px;
        s = s + alphas*ps;
        y = y + alphaz*py;
        z = z + alphaz*pz;
    end
    
    mi = sigma*mi;
end

