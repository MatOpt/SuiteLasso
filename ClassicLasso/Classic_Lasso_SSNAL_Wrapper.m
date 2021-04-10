%%*************************************************************************
%% SSNAL:
%% A semismooth Newton augmented Lagrangian method for solving Lasso problems
%% (P) min {1/2 ||Ax - b||^2 + lambda ||x||_1}
%% where lambda >0 is given data
%%
%% (D) max {-1/2 ||xi||^2  + <b,xi> | A^T xi + y = 0, ||y||_{infty} <= lambda}
%%*************************************************************************
%% SSNAL: 
%% Copyright (c) 2017 by
%% Xudong Li, Defeng Sun, and Kim-Chuan Toh
%%*************************************************************************
%% Although our code was mainly designed for the high-dimensional setting,  
%% i.e., m (#of samples) << n (#of features), we provide this wrapper file
%% to hand the case where m >> n and n <= 10,000
%%*************************************************************************

function [obj,y,xi,x,info,runhist] = Classic_Lasso_SSNAL_Wrapper(Ainput,b,n,lambda,options,y,xi,x)

isQR = 1;
QRtol = 1e-8;
if isstruct(Ainput)
   if isfield(Ainput,'A')
      A = Ainput.A; 
   else
      isQR = 0;
   end
else
   A = Ainput;
end

m = length(b);
if m > n 
   if n > 1e4
      isQR = 0;
   end
else
   isQR = 0;
end

if isQR == 0
   fprintf('\n call original SSNAL to solve the problem');
   if ~exist('x','var'); x = zeros(n,1); end;
   if ~exist('y','var'); y = zeros(n,1); end;
   if ~exist('xi','var'); xi = zeros(length(b),1); end;
   [obj,y,xi,x,info,runhist] = Classic_Lasso_SSNAL(Ainput,b,n,lambda,options,y,xi,x);
else
   fprintf('\n case m = %d, n= %d (m>n), using QR decomposition to equivalently rewirte the original problem', m,n);
   [Q,R] = qr(A,0);
   if abs(R(n,n))< QRtol
      idx = abs(diag(R)) >= QRtol;
      Q = Q(:,idx);
      R = R(idx,:);
   end
   [m,n] = size(R);
   fprintf('\n new data constructed with m = %d, n=%d', m,n);
   borg = b;
   b = (borg'*Q)'; %you may want to save new data (R,b) if solution path is needed
   const = 0.5*(norm(borg)^2  - norm(b)^2);
   options.orgojbconst = const;
   [obj,y,xi,x,info,runhist] = Classic_Lasso_SSNAL(R,b,n,lambda,options);
end
   



