%%*************************************************************************
%% SSNAL:
%% A semismooth Newton augmented Lagrangian method for solving fused Lasso problems
%% (P) min {1/2 ||Ax - b||^2 + lambda1 ||x||_1 + lambda2 ||Bx||_1}
%% where lambda1, lambda2 >0 are given data 
%% Bx = [x1 - x2;...; x(n-1) - xn];
%%*************************************************************************
%% SSNAL: 
%% Copyright (c) 2017 by
%% Xudong Li, Defeng Sun, and Kim-Chuan Toh
%%*************************************************************************
function [obj,x,xi,u,info,runhist] = Fused_Lasso_SSNAL(Ainput,b,n,lambda1,lambda2,options,x0,xi0,u0)

   rng('default');
   maxiter = 5000;
   stoptol = 1e-6;
   printyes=1;
   scale = 1;
   rescale = 1;
   precond = 2;
   stopop = 2;
   orgojbconst = 0;
   
   if isfield(options,'maxiter');  maxiter  = options.maxiter; end
   if isfield(options,'stoptol');  stoptol  = options.stoptol; end
   if isfield(options,'printyes'); printyes = options.printyes; end
   if isfield(options,'rescale'); rescale = options.rescale; end
   if isfield(options,'Lip'); Lip = options.Lip; end
   if isfield(options,'stopop'); stopop = options.stopop; end
   if isfield(options,'precond'); precond = options.precond; end
   if isfield(options,'orgojbconst'); orgojbconst = options.orgojbconst; end
%%
   fprintf('\n*************************************************************************************');
   fprintf('\n SuiteLasso');
   fprintf('\n Authors: Xudong Li, Defeng Sun, and Kim-Chuan Toh                                ');
   fprintf('\n*************************************************************************************');

   %% Amap and ATmap
%%
   tstart = clock;
   m = length(b);
   if isstruct(Ainput)
      if isfield(Ainput,'A')
         A = Ainput.A; 
         Amap0 = @(x) A*x;      
      elseif isfield(Ainput,'Amap')
         Amap0 = Ainput.Amap; 
      else
         error(' Ainput.Amap not defined')         
      end
      if exist('A','var')
         ATmap0 = @(y) A'*y;       
      elseif isfield(Ainput,'ATmap')
         ATmap0 = Ainput.ATmap;    
      else
         error(' Ainput.ATmap not defined')
      end
   else
      A = Ainput; 
      Amap0 = @(x) A*x;
      ATmap0 = @(y) A'*y;
   end
   if ~exist('A','var')
      fprintf('\n Remark: better performance can be achieved if the'); 
      fprintf('\n matrix representation of the linear operator A is available');
   end
   
   AATmap0 = @(x) Amap0(ATmap0(x));
   if ~exist('Lip','var')
      eigsopt.issym = 1;
      tstartLip = clock;
      Lip = eigs(AATmap0,m,1,'LA',eigsopt);
      fprintf('\n Lip = %3.2e, time = %3.2f', Lip, etime(clock, tstartLip));
   end
   
   [Bmap0,BTmap0] = FLBmap(n);
   normb = 1 + norm(b);
   if ~exist('x0','var') || ~exist('xi0','var') || ~exist('u0','var') 
      x = zeros(n,1); xi = zeros(m,1); u = zeros(n,1);
   else
      x = x0; xi = xi0; u = u0;
   end
   
   %%
   parmain.m = m; parmain.n = n;
   parmain.scale = scale;
   parmain.tstart = tstart;
   parmain.rescale = rescale;
   parmain.stoptol = stoptol;
   parmain.Lip = Lip;
   if isfield(options,'sigma'); parmain.sigma = options.sigma; end
   parmain.maxiter = maxiter;
   parmain.stopop = stopop;
   parmain.precond = precond;
   parmain.printyes = printyes;
   parmain.orgojbconst = orgojbconst;
   if exist('A','var'); parmain.A = A; end
   
   [obj_main,u,xi,x,info_main,runhist_main] = ...
        Fused_Lasso_SSNAL_main(Amap0,ATmap0,Bmap0,BTmap0,b,lambda1,lambda2,parmain,u,xi,x);
    
   iter = info_main.iter;
   bscale = info_main.bscale;
   cscale = info_main.cscale;
   Ax = info_main.Ax;
   ttime = info_main.ttime;
   msg = info_main.msg;
   %% 
%%-----------------------------------------------------------------
%% recover orignal variables
%%-----------------------------------------------------------------
   if (iter == maxiter)
      msg = ' maximum iteration reached';
      info.termcode = 3;
   end   
   xi = xi*sqrt(bscale*cscale);
   Atxi = ATmap0(xi);
   u = u*cscale;
   x = x*bscale;
   %up = up*bscale;
   Ax = Ax*sqrt(bscale*cscale);
   Bx = Bmap0(x);
   Rd = Atxi + u;
   Rp1 = Ax - b;
   normRp = norm(Rp1 - xi);
   normRp1 = norm(Rp1);
   normRd = norm(Rd);
   normu = norm(u);
   primfeasorg = normRp/normb;
   dualfeasorg = normRd/(1+normu);
   primobj = 0.5*normRp1^2 + lambda1*norm(x,1) + lambda2*norm(Bx,1) + orgojbconst;
   dualobj = -(0.5*norm(xi)^2 + b'*xi) + orgojbconst;
   relgap = (primobj-dualobj)/( 1+abs(primobj)+abs(dualobj));
   obj = [primobj, dualobj];
   grad = ATmap0(Rp1);
   etaorg = errcom(x,grad,lambda1,lambda2);
   %etaorg = norm(x- proxFL(Binput,x - grad,lambda1org,lambda2org));
   eta = etaorg/(1+norm(grad)+norm(x));
   runhist.m = m;
   runhist.n = n;
   
   runhist.iter = iter;
   runhist.totaltime = ttime;
   runhist.primobjorg = primobj; 
   runhist.dualobjorg = dualobj;
   runhist.maxfeas = max([dualfeasorg, primfeasorg]);
   runhist.eta = eta;
   runhist.etaorg = etaorg;
   info.m = m;
   info.n = n;
   info.minx = min(min(x));
   info.maxx = max(max(x));
   info.relgap = relgap;
   info.iter = iter;
   info.time = ttime;
   info.eta = eta;
   info.etaorg = etaorg;
   info.obj = obj;
   info.dualfeasorg = dualfeasorg;
   info.primfeasorg = primfeasorg;
   info.maxfeas = max([dualfeasorg, primfeasorg]);
   info.Axmb = normRp1;
   info.nnzx = findnnz(x,0.999);
   info.nnzBx = findnnz(Bx,0.999);
   if (printyes) 
      if ~isempty(msg); fprintf('\n  %s',msg); end
      fprintf('\n--------------------------------------------------------------');
      fprintf('------------------');
      fprintf('\n  number iter = %2.0d',iter);       
      fprintf('\n  time = %3.2f',ttime);       
      fprintf('\n  time per iter = %5.4f',ttime/iter); 
      fprintf('\n  primobj = %10.9e, dualobj = %10.9e, relgap = %3.2e',primobj,dualobj, relgap);    
      fprintf('\n  primfeasorg    = %3.2e, dualfeasorg    = %3.2e',...
	      primfeasorg, dualfeasorg); 
      fprintf('\n  eta = %3.2e, etaorg = %3.2e', eta, etaorg);
      fprintf('\n  min(X)    = %3.2e, max(X)    = %3.2e',...
          info.minx,info.maxx); 
      fprintf('\n  number of nonzeros in x (0.999)  = %3.0d', findnnz(x,0.999));
      fprintf('\n  number of nonzeros in Bx (0.999) = %3.0d', findnnz(Bx,0.999));
      fprintf('\n--------------------------------------------------------------');
      fprintf('------------------\n');
   end
%%**********************************************************************


   