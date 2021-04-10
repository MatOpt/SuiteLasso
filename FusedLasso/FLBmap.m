function [Bmap,BTmap] = FLBmap(n)

%B = spdiags([ones(n,1),-ones(n,1)],[0,1]',n-1,n);
%BBt = spdiags([-ones(n-1,1),2*ones(n-1,1),-ones(n-1,1)],[-1,0,1]',n-1,n-1);
Bmap = @(x) x(1:end-1) - x(2:end);
BTmap = @(y) [y;0] - [0;y];