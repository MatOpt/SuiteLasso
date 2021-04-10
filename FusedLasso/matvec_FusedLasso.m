

function FLxi = matvec_FusedLasso(xi,par,~)

AP = par.AP;
PU1 = par.PU1;
if (size(PU1,2) > 0)
   tmp = (xi'*AP)';  
   tmp2 = (tmp'*PU1)';
   FLxi = xi + par.sigma*(AP*(par.Ph.*tmp + PU1*tmp2));
elseif (norm(par.Ph) > 0)
   tmp = (xi'*AP)';      
   FLxi = xi + par.sigma*(AP*(par.Ph.*tmp));
else
   FLxi = xi; 
end