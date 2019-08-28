%xini = obj.xall;
%fun = @(x)(fufu);


function [ret] = fufu(obj,x)
obj.xall = x;
obj.eval_alx();
ret = obj.ALx;
end