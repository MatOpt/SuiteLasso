%%*********************************************************************
%%*********************************************************************
clear all;
rng('default');
addpath(genpath(pwd));

test_type = 1;
%% Input data: A b lambda 
if test_type == 0
 m = 50;
 n = 50000;
 A = randn(m,n);
 k = 0.1*n;
 xx = zeros(n,1);
 xstar = 1.2*sqrt(2*log(n));
 xx(1:k) = xstar;
 %err = 0.001*randn(m,1);
 b = A*xx;%+err;
elseif test_type == 1
   m = 50000;
   n = 50;
   A = rand(m,n);
   b = rand(m,1);
elseif test_type == 2
   Datapath = [pwd,filesep,'UCIdata',filesep,'uci_CT.mat'];
   data = load(Datapath);
   A = data.A;
   b = data.b;
   [m,n] = size(A);
end
 Amap  = @(x) A*x;
 ATmap = @(x) A'*x; 
    
%% tuning parameters   

fprintf('\n-----------------------------------------------');
fprintf('------------------------------')
fprintf('\n Problem: n = %g,  m = %g  ',n,m)
fprintf('\n-----------------------------------------------');
fprintf('------------------------------')

stoptol = 1e-6;%stopping tol

lambdaop = 2;
if lambdaop == 1
   lambda = 1e-3;
else
   lambdamax=norm(A'*b,'inf');
   fprintf('\n lambdamax = %3.2e', lambdamax);
   lambda=1e-3*lambdamax;
end

SSNAL = 1;
if SSNAL
   nalop.stoptol = stoptol;
   nalop.Lip = 1;
   Ainput.A = A;
   [obj,y,xi,x,info,runhist] = Classic_Lasso_SSNAL_Wrapper(Ainput,b,n,lambda,nalop);
end