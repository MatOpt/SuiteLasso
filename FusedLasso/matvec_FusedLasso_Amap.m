

function FLxi = matvec_FusedLasso_Amap(xi,par,Ainput)

PU1 = par.PU1;
tmp = Ainput.ATmap(xi);
tmp = tmp(par.info_u.rr1);

FLxi = xi;
if (size(PU1,2) >0) || (norm(par.Ph) > 0)
   tmp4 = zeros(par.n,1);
   if norm(par.Ph) > 0
      tmp3 = par.Ph.*tmp;
   end
   if size(PU1,2) > 0
      tmp2 = (tmp'*PU1)';
      if exist('tmp3','var') 
         tmp3 = tmp3 + PU1*tmp2;
      else
         tmp3 = PU1*tmp2;
      end
   end
   tmp4(par.info_u.rr1) = tmp3;
   FLxi = FLxi + par.sigma*Ainput.Amap(tmp4);
end
% if (size(PU1,2) > 0)
%    tmp2 = (tmp'*PU1)';
%    tmp3 = par.Ph.*tmp + PU1*tmp2;
%    tmp4 = zeros(par.n,1);
%    tmp4(par.info_u.rr1) = tmp3;
%    FLxi = xi + par.sigma*(par.Amap(tmp4));
% elseif (norm(par.Ph) > 0)
%    tmp3 = par.Ph.*tmp;
%    tmp4 = zeros(par.n,1);
%    tmp4(par.info_u.rr1) = tmp3;
%    FLxi = xi + par.sigma*(par.Amap(tmp4));
% else
%    FLxi = xi; 
% end