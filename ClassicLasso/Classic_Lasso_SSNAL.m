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
function [obj,y,xi,x,info,runhist] = Classic_Lasso_SSNAL(Ainput,b,n,lambda,options,y,xi,x)

   rng('default');
   maxiter = 5000;
   stoptol = 1e-6;
   printyes = 1;
   scale = 0;
   dscale = ones(n,1);
   rescale = 1;
   Lip = 1;
   Ascale = 0;
   Ascaleyes = 0;
   orgojbconst = 0;
   
   if isfield(options,'maxiter');  maxiter  = options.maxiter; end
   if isfield(options,'stoptol');  stoptol  = options.stoptol; end
   if isfield(options,'printyes'); printyes = options.printyes; end
   if isfield(options,'rescale'); rescale = options.rescale; end
   if isfield(options,'Lip'); Lip = options.Lip; end
   if isfield(options,'Ascale'); Ascale = options.Ascale; end
   if isfield(options,'orgojbconst'); orgojbconst = options.orgojbconst; end
%%
   fprintf('\n*************************************************************************************');
   fprintf('\n SuiteLasso');
   fprintf('\n Authors: Xudong Li, Defeng Sun, and Kim-Chuan Toh');
   fprintf('\n*************************************************************************************');
   tstart = clock;
   m = length(b);
   if isstruct(Ainput)
      if isfield(Ainput,'A'); A = Ainput.A; end
      if isfield(Ainput,'Amap'); Amap0 = Ainput.Amap; end
      if isfield(Ainput,'ATmap'); ATmap0 = Ainput.ATmap; end
   else
      A = Ainput; 
      Amap0 =@(x) A*x;
      ATmap0 = @(y) A'*y;
   end
   existA = exist('A','var');
   if Ascale ~= 0 && existA
      tscale = clock; 
      if Ascale == 1
         dscale = 1./max(1,sqrt(sum(A.*A))');
      elseif Ascale == 2
         dscale = 1./sqrt(sum(A.*A))';
      end
      A = A*spdiags(dscale,0,n,n);
      Ascaleyes = 1;
      fprintf('\n time for scaling A  = %3.1f', etime(clock, tscale));
   end      
   if ~exist('x','var') || ~exist('xi','var') || ~exist('y','var')
      x = zeros(n,1); xi = zeros(m,1); y = x;
   end
   parmain.dscale    = dscale;
   parmain.Ascaleyes = Ascaleyes;
   parmain.m = m; parmain.n = n;
   parmain.scale = scale;
   parmain.existA = existA;
   parmain.orgojbconst = orgojbconst;
   
   if existA; parmain.A = A; end
   %parmain.Amap0 = Amap0; parmain.ATmap0 = ATmap0;
   
   parmain.tstart = tstart;   
   parmain.Lip = Lip;
   if isfield(options,'sigma'); parmain.sigma = options.sigma; end
   parmain.maxiter = maxiter;
   parmain.printyes = printyes;
   parmain.rescale = rescale;
   parmain.stoptol = stoptol;
   
   [obj_main,y,xi,x,info_main,runhist] = ...
       Classic_Lasso_SSNAL_main(Amap0,ATmap0,b,lambda,parmain,y,xi,x);   
   iter   = info_main.iter;
   bscale = info_main.bscale;
   cscale = info_main.cscale;
   Ax = info_main.Ax;
   ttime = info_main.ttime;
   msg = info_main.msg;
%%-----------------------------------------------------------------
%% recover orignal variables
%%-----------------------------------------------------------------
   if (iter == maxiter)
      msg = ' maximum iteration reached';
   end
   xi = xi*sqrt(bscale*cscale);
   Atxi = ATmap0(xi);
   y = y*cscale;
   x = x*bscale;
   if Ascaleyes; x = dscale.*x; y = y./dscale; end
   Rd = Atxi + y;
   dualfeasorg = norm(Rd)/(1+norm(y));
   Ax = Ax*sqrt(bscale*cscale);
   Rp = Ax - b + xi;
   primfeasorg = norm(Rp)/(1 + norm(b));
   primobj = 0.5*norm(Ax - b)^2 + lambda*norm(x,1) + orgojbconst;
   dualobj = -0.5*norm(xi)^2 + b'*xi + orgojbconst;
   relgap = (primobj-dualobj)/( 1+abs(primobj)+abs(dualobj));
   obj = [primobj, dualobj];
   grad = ATmap0(Ax - b);
   etaorg = norm(grad + proj_inf(x - grad,lambda));  
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
   info.maxfeas = max([dualfeasorg, primfeasorg]);
   info.cnt_Amap = sum(runhist.cnt_Amap);
   info.cnt_ATmap = sum(runhist.cnt_ATmap);
   info.cnt_pAATmap = sum(runhist.cnt_pAATmap);
   info.cnt_fAATmap = sum(runhist.cnt_fAATmap);
   info.nnz = findnnz(x,0.999);
   info.x  = x;
   if (printyes) 
      if ~isempty(msg); fprintf('\n  %s',msg); end
      fprintf('\n--------------------------------------------------------------');
      fprintf('------------------');
      fprintf('\n  number iter = %2.0d',iter);       
      fprintf('\n  time = %3.2f',ttime);       
      fprintf('\n  time per iter = %5.4f',ttime/iter); 
      fprintf('\n  primobj = %9.8e, dualobj = %9.8e, relgap = %3.2e',primobj,dualobj, relgap); 
      fprintf('\n  primfeasorg    = %3.2e, dualfeasorg    = %3.2e',...
	      primfeasorg, dualfeasorg); 
      fprintf('\n  eta = %3.2e, etaorg = %3.2e', eta, etaorg);
      fprintf('\n  min(X) = %3.2e, max(X) = %3.2e',...
          info.minx,info.maxx); 
      fprintf('\n  Amap cnt = %3d, ATmap cnt = %3d, partial AATmap cnt = %3d, full AATmap cnt = %3d',...
          info.cnt_Amap, info.cnt_ATmap, info.cnt_pAATmap, info.cnt_fAATmap);
      fprintf('\n  number of nonzeros in x (0.999) = %3.0d',findnnz(x,0.999)); 
      fprintf('\n--------------------------------------------------------------');
      fprintf('------------------\n');
   end
%%**********************************************************************

 
   