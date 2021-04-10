function [k,xnew] = cardcal(x,r,tiny)

if ~exist('tiny','var'); tiny = 1e-16; end

n = length(x);
normx1 = norm(x,1);

if min([normx1,norm(x,'inf')]) <= tiny
   k = 0;
   if nargout >1
      xnew = zeros(n,1);
   end
   return
end

[absx,idx] = sort(abs(x),'descend');
tmpidx = find(cumsum(absx) > r*normx1);
k = tmpidx(1);

if nargout > 1
   xnew = zeros(n,1);
   idxnew = idx(1:k);
   xnew(idxnew) = x(idxnew);
end