function [E] = E(x,s,y,z,mi,p,fx,hx,gx,h,g)

a = norm( fx(x) - hx(x)'*y - gx(x)'*z ,p);
b = norm( diag(s)*z - mi*ones(length(s),1) ,p);
c = norm( h(x) ,p);
d = norm( g(x) - s ,p);

E = max([a,b,c,d])

end

