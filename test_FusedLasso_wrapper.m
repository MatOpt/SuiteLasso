%%*********************************************************************
%%*********************************************************************
clear all;
rng('default');
addpath(genpath(pwd));
 
%% Input data: A b lambda1 lambda2 
testtype = 1; 

if testtype == 1
   n = 50;
   m = 50000;
else
   n = 50000;
   m = 50;
end
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

fprintf('\n-----------------------------------------------');
fprintf('------------------------------')
fprintf('\n Problem: n = %g,  m = %g    lambda(max) = %g ',n,m, lambdamax)
fprintf('\n-----------------------------------------------');
fprintf('------------------------------')

stoptol = 1e-6;%stopping tol

for crho = 1e-3 
   lambda1 = crho*lambdamax;
   lambda2 = 2*lambda1;
   if (true)
      opts.stoptol = stoptol;
      opts.Lip = 1; %can be set to 1
      Ainput.A = A;
      Ainput.Amap = @(x) Amap(x);
      Ainput.ATmap = @(x) ATmap(x);
      [obj,x,xi,u,info,runhist] = ...
              Fused_Lasso_SSNAL_Wrapper(Ainput,b,n,lambda1,lambda2,opts);
   end
end
%%*********************************************************************
