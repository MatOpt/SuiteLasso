function Ay = matvec_ClassicLasso(y,par,AP)
tmp = AP'*y;
Ay = y + par.sigma*(AP*tmp);