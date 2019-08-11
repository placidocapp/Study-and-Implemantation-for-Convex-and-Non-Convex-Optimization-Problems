function [fi] = Dmerit(x,s,mi,v,f,fx,h,hx,g,gx,px,ps,p)
% Directional derivative of merit function with respect to (px,ps)

if p == 1
    fi = fx(x)'*px - mi*sum(log(s)'*ps) + v*( hx(x)'*sign(h(x)) )'*px + ...
         v*( gx(x)'*sign(g(x)) )'*px -  ...
         v*( ones(length(s),1)*norm(g(x)-s,p) )'*ps;
elseif p == 2
    fi = fx(x)'*px - mi*sum(log(s)'*ps) + v*( 2*hx(x)'*h(x) )*px + ...
         v*( 2*gx(x)'*g(x) )*px -  ...
         v*( ones(length(s),1)'*norm(g(x)-s,p) )*ps;
else
    disp('p não é valido, escolha 1 ou 2')

end