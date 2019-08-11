function [fi] = merit(x,s,mi,v,f,h,g,p)
% Merit function

fi = f(x) - mi*sum(log(s)) + v*norm(h(x),p) + v*norm(g(x) - s,p);

end


