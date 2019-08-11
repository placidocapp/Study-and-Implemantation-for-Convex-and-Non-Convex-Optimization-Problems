function [val] = mul(z,H)
val = 0;
for i = 1:length(z)
   val = val + z(i)*H(:,:,i);
end
end

