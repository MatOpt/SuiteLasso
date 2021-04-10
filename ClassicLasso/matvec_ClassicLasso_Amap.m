function Ay = matvecIpxi(y,par,Ainput)
tmp = Ainput.ATmap(y);
tmp = (1 - par.rr).*tmp;
Ay = y + par.sigma*Ainput.Amap(tmp);