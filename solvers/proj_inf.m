function [y,rr] = proj_inf(x,lambda)
if lambda < 0; 
   error('Lambda needs to be nonegative');
elseif lambda == 0; 
   y = 0*x; 
else
   y = max(-lambda,min(x,lambda));
end
rr = (y == x); %(abs(y-x) < 1e-6*lambda);