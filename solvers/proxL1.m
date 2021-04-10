%% proximal mapping for L1 norm 
%% y = argmin{ 0.5*||y - x||^2 + \lmabda ||y||_1}
%%
function [y,rr] = proxL1new(x,lambda)
[tmp, rr] = proj_inf(x,lambda);
y = x - tmp;
rr = ~rr;
%y = x - proj_inf(x,lambda);