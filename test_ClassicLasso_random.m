%%*********************************************************************
%%*********************************************************************
clear all;
rng('default');
addpath(genpath(pwd));
 
%% Input data: A b lambda 
 
 m = 50;
 n = 50000;
 A = randn(m,n);
 k = 0.1*n;
 xx = zeros(n,1);
 xstar = 1.2*sqrt(2*log(n));
 xx(1:k) = xstar;
 %err = 0.001*randn(m,1);
 b = A*xx;%+err;
 
 Amap  = @(x) A*x;
 ATmap = @(x) A'*x; 
    
%% tuning parameters   

lambdamax=norm(ATmap(b),'inf');
eigsopt.issym = 1;
Rmap=@(x) A*x;
Rtmap=@(x) A'*x;
RRtmap=@(x) Rmap(Rtmap(x));
Lip = eigs(RRtmap,length(b),1,'LA',eigsopt);

fprintf('\n-----------------------------------------------');
fprintf('------------------------------')
fprintf('\n Problem: n = %g,  m = %g    lambda(max) = %g ',n,m, lambdamax)
fprintf('\n Lip = %3.2e', Lip);
fprintf('\n-----------------------------------------------');
fprintf('------------------------------')

stoptol = 1e-6;%stopping tol

for crho = 1e-3 
   lambda=crho*lambdamax;
   if (true)
      opts.stoptol = stoptol;
      opts.Lip = Lip; %can be set to 1
      Ainput.A = A;
      Ainput.Amap = @(x) Amap(x);
      Ainput.ATmap = @(x) ATmap(x);
      if (false) %% supply initial point
         x0 = zeros(n,1);
         y0 = zeros(n,1);
         xi0= zeros(m,1);      
         [obj,y,xi,x,info,runhist] = Classic_Lasso_SSNAL(Ainput,b,n,lambda,nalop,y0,xi0,x0);
      else %% no initial point available
         [obj,y,xi,x,info,runhist] = Classic_Lasso_SSNAL(Ainput,b,n,lambda,opts);                                     
      end
   else
      [obj,y,xi,x,info,runhist] = Classic_Lasso_ADMM(Ainput,b,n,lambda,opts);
   end
end
%%*********************************************************************
