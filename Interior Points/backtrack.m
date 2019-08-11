function [alpha] = backtrack(x,s,px,ps,ro,po,p,alpha,f,fx,h,hx,g,gx,mi,c1)

v = (fx(x)'*px + mi*sum(log(s))) / ((1-po)*(norm(h(x),p)+norm(g(x)-s,p)));

%Reduce the step until it satisfies the Armijo condition 
while merit(x,s,mi,v,f,h,g,p) > merit(x,s,mi,v,f,h,g,p) + ...
        c1*alpha*g(x)'*Dmerit(x,s,mi,v,f,fx,h,hx,g,gx,px,ps,p) 
    alpha = ro*alpha;
end

end

