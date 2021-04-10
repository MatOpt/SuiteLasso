%%*************************************************************************
%% psqmr:  preconditioned symmetric QMR with left (symmetric) preconditioner. 
%%
%% b = rhs vector.
%% resnrm = norm of qmr-generated residual vector b-Ax. 
%%
%% SDPNAL: 
%% Copyright (c) 2008 by
%% Xinyuan Zhao, Defeng Sun, and Kim-Chuan Toh 
%%*************************************************************************

   function  [x,Ax,resnrm,solve_ok] = psqmry(matvecfname,A,b,par,x0,Ax0) 

   N = length(b); 
   maxit = max(5000,sqrt(length(b))); 
   tol = 1e-6*norm(b);
   stagnate_check = 20;
   miniter = 0;
   
   if (nargin < 5); x0 = zeros(N,1); end
   if isfield(par,'maxit'); maxit = par.maxit; end
   if isfield(par,'tol');  tol = par.tol; end
   if isfield(par,'stagnate_check_psqmr') 
      stagnate_check = par.stagnate_check_psqmr; 
   end
   if isfield(par,'minitpsqmr'); miniter = par.minitpsqmr; end

   solve_ok = 1; 
   printlevel = 0; 
%%
   x = x0; 
   if (norm(x) > 0)
      if nargin < 6; Ax0 = feval(matvecfname,x0,par,A); end
      Aq = Ax0;
   else
      Aq = zeros(N,1);  
   end
   r = b-Aq;  
   err = norm(r); resnrm(1) = err; minres = err; 
   %if err < tol; Ax = Ax0; return; end
%%
   q = precondfun(par,r); 
   tau_old  = norm(q);      
   rho_old  = r'*q; 
   theta_old = 0; 
   d = zeros(N,1); 
   res = r; Ad = zeros(N,1);
%%      
%% main loop
%%
   tiny = 1e-30; 
   for iter = 1:maxit 
       Aq = feval(matvecfname,q,par,A);
       sigma = q'*Aq; 
       if (abs(sigma) < tiny)
          solve_ok = 2; 
          if (printlevel); fprintf('s1'); end
          break;
       else
          alpha = rho_old/sigma; 
          r = r - alpha*Aq;
       end
       u = precondfun(par,r); 
       %%
       theta = norm(u)/tau_old; c = 1/sqrt(1+theta^2); 
       tau = tau_old*theta*c;
       gam = (c^2*theta_old^2); eta = (c^2*alpha); 
       d = gam*d + eta*q;
       x = x + d; 
       %%----- stopping conditions ----
       Ad = gam*Ad + eta*Aq;
       res = res - Ad; 
       err = norm(res); resnrm(iter+1) = err; 
       if (err < minres); minres = err; end
       if (err < tol) & (iter > miniter) & (b'*x > 0); break; end  
       if (iter > stagnate_check) & (iter > 10)
          ratio = resnrm(iter-9:iter+1)./resnrm(iter-10:iter); 
          if (min(ratio) > 0.997) & (max(ratio) < 1.003)
             if (printlevel); fprintf('s'); end
             solve_ok = -1; 
             break;
          end       
       end
       %%----------------------------- 
       if (abs(rho_old) < tiny)
          solve_ok = 2; 
          fprintf('s2');
          break;
       else
          rho  = r'*u; 
          beta = rho/rho_old; 
          q = u + beta*q; 
       end
       rho_old = rho; 
       tau_old = tau; 
       theta_old = theta; 
   end
   if (iter == maxit); solve_ok = -2; end
   if (solve_ok ~= -1); 
      if (printlevel); fprintf(' '); end
   end
   Ax = b -res;
%%*************************************************************************

   function  q = precondfun(par,r) 

   precond = 0; 
   if isfield(par,'precond'); precond = par.precond; end
   if (precond == 0)
      q = r; 
   elseif (precond == 1)
      q = par.invdiagM.*r;
   elseif (precond == 2)
       dlast = par.d(end);
       tmp1 = 2/dlast; 
       tmp2 = 1./par.d - tmp1; 
       q = tmp1*r + par.V*(tmp2.*(par.Vt*r));
   elseif (precond == 3)
      q = par.invM(r);
   elseif (precond==4)
      q = mylinsysolve(par.L,r);   
   end
%%*************************************************************************
