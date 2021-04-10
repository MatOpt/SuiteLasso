function [u,Atxi,xi,par,runhist,info] = ...
    Fused_Lasso_SSNCG(b,Ainput,Binput,x0,Ax0,Atxi0,xi0,...
    lambda1,lambda2,par,options)

printsub = 1;
breakyes = 0;
maxitersub = 50;
tiny = 1e-10;
tol = 1e-6;
maxitpsqmr =500;
use_proximal = 0;

if isfield(options,'printsub'); printsub = options.printsub; end
if isfield(options,'maxitersub'); maxitersub = options.maxitersub; end
if isfield(options,'tiny'); tiny = options.tiny; end
if isfield(options,'tol'); tol = min(tol,options.tol); end
if isfield(options,'use_proximal'); use_proximal = options.use_proximal; end

sig = par.sigma;
bscale = options.bscale;
cscale = options.cscale;
normborg = 1+norm(b)*sqrt(bscale*cscale);
%% preperation
Amap = @(x) Ainput.Amap(x);
ATmap = @(x) Ainput.ATmap(x);
%par.lsAmap = Amap;
%Bmap = Binput.Bmap;
%BTmap = Binput.BTmap;
uinput = x0 - sig*Atxi0;
[up,info_u] = proxFL(Binput,uinput,sig*lambda1,sig*lambda2);
par.info_u = info_u;
par.innerNT = par.innerNT + info_u.innerNT;
par.innerflsa = par.innerflsa + info_u.innerflsa; 
u = (uinput -up)/sig;
Rp = Ax0 - b - xi0;
normRp = norm(Rp);
Atxi = Atxi0; xi = xi0;
Ly = -b'*xi - 0.5*norm(xi)^2  - norm(up)^2/(2*sig);
%par.cginitialxi = xi;
runhist.psqmr(1) = 0;
runhist.findstep(1) = 0;
%% main Newton iteration
for itersub = 1:maxitersub    
    xiold = xi; Atxiold = Atxi;
    Rd = Atxi + u;
    normRd = norm(Rd); 
    %Ly = -b'*xi - 0.5*norm(xi)^2  - norm(up)^2/(2*sig);
    Aup = Amap(up);
    GradLxi = -(xi + b - Aup);
    normGradLxi = norm(GradLxi)*sqrt(bscale*cscale)/normborg;
    priminf_sub = normGradLxi; 
    dualinf_sub = normRd*cscale/(1+norm(u)*cscale);
    if max(priminf_sub,dualinf_sub) < tol
       tolsubconst = 0.1;
    else
       tolsubconst = 0.05; %1e-2;
    end
    tolsub = max(min(1,par.tolconst*dualinf_sub),tolsubconst*tol);
    runhist.priminf(itersub) = priminf_sub;
    runhist.dualinf(itersub) = dualinf_sub;
    runhist.Ly(itersub)      = Ly;
    if (printsub)
        fprintf('\n      %2.0d  %- 11.10e %3.2e %3.2e %3.2e',...
                 itersub,Ly,priminf_sub,dualinf_sub,par.tolconst);
    end
    if max(normGradLxi) < max(tolsub) && itersub > 1
        msg = 'good termination in subproblem:';
        if printsub
            fprintf('\n       %s  ',msg);
            fprintf(' normRd=%3.2e, gradLyxi = %3.2e, tolsub=%3.2e',...
                          normRd,normGradLxi,tolsub);
        end
        breakyes = -1;
        break;
%     elseif (max(priminf_sub,dualinf_sub) < 0.5*tol)
%        msg = sprintf('max(priminf_sub,dualinf_sub) < %3.2e',0.5*tol);
%        fprintf('\n       %s',msg);
%        breakyes = -1;
%        break;
    end
    %% Compute Newton direction
    %% precond = 0, 
     par.epsilon = min([1e-3,0.1*normGradLxi]); %% good to add
     if (dualinf_sub > 1e-3) || (itersub <= 5)
          maxitpsqmr = max(maxitpsqmr,200); 
     elseif (dualinf_sub > 1e-4)	 
          maxitpsqmr = max(maxitpsqmr,300); 
     elseif (dualinf_sub > 1e-5)	 
          maxitpsqmr = max(maxitpsqmr,400); 
     elseif (dualinf_sub > 5e-6)
          maxitpsqmr = max(maxitpsqmr,500); 
     end
     if (itersub > 1) 
          prim_ratio = priminf_sub/runhist.priminf(itersub-1); 
          dual_ratio = dualinf_sub/runhist.dualinf(itersub-1); 
     else
          prim_ratio = 0; dual_ratio = 0;
     end
     rhs = GradLxi;
     tolpsqmr = min([5e-3, 0.1*norm(rhs)]);
     const2 = 1;
     if itersub > 1 && (prim_ratio > 0.5 || priminf_sub > 0.1*runhist.priminf(1))
         const2 = 0.5*const2;
     end
     if (dual_ratio > 1.1); const2 = 0.5*const2; end
     tolpsqmr = const2*tolpsqmr;
     par.tol = tolpsqmr; 
     par.maxit = 2*maxitpsqmr;
     [dxi,resnrm,solve_ok,par] = Fused_Lasso_linsys_solver(Ainput,rhs,par);
     Atdxi = ATmap(dxi);
     iterpsqmr = length(resnrm)-1;
     if (printsub)
        fprintf('| %3.1e %3.1e %3.0d',par.tol,resnrm(end),iterpsqmr);
        fprintf(' %2.1f [%2d, %2d, (%2d, %2d)]',const2,sum(par.info_u.rr1),sum(1-par.info_u.rr2), par.lenP, par.numblk1);
     end
     par.iter = itersub;
     if (itersub<=3) && (dualinf_sub > 1e-4) || (par.iter <3)
         stepop = 1;
     else
         stepop = 2;
     end
     steptol = 1e-5; 
     step_op.stepop=stepop;
     [par,Ly,xi,Atxi,u,up,alp,iterstep] = ...
         findstep(par,Binput,b,lambda1,lambda2,Ly,xi,Atxi,...
            u,up,dxi,Atdxi,steptol,step_op); 
     runhist.solve_ok(itersub) = solve_ok;
     runhist.psqmr(itersub)    = iterpsqmr; 
     runhist.findstep(itersub) = iterstep; 
     Ly_ratio = 1; 
     if (itersub > 1)
          Ly_ratio = (Ly-runhist.Ly(itersub-1))/(abs(Ly)+eps);
     end
     if (printsub)
        fprintf(' %3.2e %2.0f',alp,iterstep);
        if (Ly_ratio < 0); fprintf('-'); end
     end
     %% check for stagnation
     if (itersub > 4)
         idx = [max(1,itersub-3):itersub]; 
         tmp = runhist.priminf(idx); 
         ratio = min(tmp)/max(tmp);
         if (all(runhist.solve_ok(idx) <= -1)) && (ratio > 0.9) ... 
              && (min(runhist.psqmr(idx)) == max(runhist.psqmr(idx))) ...
              && (max(tmp) < 5*tol)
              fprintf('#')
              breakyes = 1; 
         end
         const3 = 0.7;
         priminf_1half  = min(runhist.priminf(1:ceil(itersub*const3))); 
         priminf_2half  = min(runhist.priminf(ceil(itersub*const3)+1:itersub));
         priminf_best   = min(runhist.priminf(1:itersub-1)); 
         priminf_ratio  = runhist.priminf(itersub)/runhist.priminf(itersub-1);
         dualinf_ratio  = runhist.dualinf(itersub)/runhist.dualinf(itersub-1);  
         stagnate_idx   = find(runhist.solve_ok(1:itersub) <= -1);
         stagnate_count = length(stagnate_idx); 
         idx2 = [max(1,itersub-7):itersub]; 
         if (itersub >= 10) && all(runhist.solve_ok(idx2) == -1) ... 
              && (priminf_best < 1e-2) && (dualinf_sub < 1e-3)                    
              tmp = runhist.priminf(idx2); 
              ratio = min(tmp)/max(tmp); 
              if (ratio > 0.5) 
                  if (printsub); fprintf('##'); end
                  breakyes = 2; 
              end
         end
          if (itersub >= 15) && (priminf_1half < min(2e-3,priminf_2half)) ...
               && (dualinf_sub < 0.8*runhist.dualinf(1)) && (dualinf_sub < 1e-3) ...
               && (stagnate_count >= 3) 
               if (printsub); fprintf('###'); end
               breakyes = 3;
          end
          if (itersub >= 15) && (priminf_ratio < 0.1) ...
               && (priminf_sub < 0.8*priminf_1half) ...
               && (dualinf_sub < min(1e-3,2*priminf_sub)) ...
               && ((priminf_sub < 2e-3) || (dualinf_sub < 1e-5 && priminf_sub < 5e-3)) ...
               && (stagnate_count >= 3) 
               if (printsub); fprintf(' $$'); end
               breakyes = 4;  
          end
          if (itersub >=10) && (dualinf_sub > 5*min(runhist.dualinf)) ...
               && (priminf_sub > 2*min(runhist.priminf)) %% add: 08-Apr-2008
               if (printsub); fprintf('$$$'); end
               breakyes = 5;
          end
          if (itersub >= 20)
             %% add: 12-May-2010
             dualinf_ratioall = runhist.dualinf(2:itersub)./runhist.dualinf(1:itersub-1);
             idx = find(dualinf_ratioall > 1); 
             if (length(idx) >= 3)
                  dualinf_increment = mean(dualinf_ratioall(idx)); 
                  if (dualinf_increment > 1.25)
                      if (printsub); fprintf('^^'); end
                      breakyes = 6;                            
                  end                    
             end              
          end 
          if (breakyes > 0)
             Rd =  Atxi + u;
             normRd = norm(Rd);
             Aup = Amap(up);
             fprintf('\n new dualfeasorg = %3.2e', normRd*cscale/(1+norm(u)*cscale));
             break
          end
     end 
end
info.breakyes = breakyes;
info.itersub = itersub;
info.tolconst = par.tolconst;
info.up = up;
info.Aup = Aup;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [par,xi,Atxi,u,up,alp,iter] = ...
%          findstep(par,Binput,b,x0,lambda1,lambda2,xi0,Atxi0,...
%          u0,up0,dxi,Atdxi,tol,options)
% 
%    printlevel = 0; 
%    maxit = max(20,ceil(log(1/(tol+eps))/log(2)));
%    c1 = 1e-4;   
%    sig = par.sigma;
% %%  
%    g0  = dxi'*(-b - xi0) + Atdxi'*up0;
%    Ly0 = -b'*xi0 - 0.5*norm(xi0)^2 - norm(up0)^2/(2*sig);
%    if (g0 <= 0)
%       alp = 0; iter = 0; 
%       if (printlevel); 
%          fprintf('\n Need an ascent direction, %2.1e  ',g0); 
%       end
%       xi = xi0;
%       Atxi = Atxi0;
%       u = u0;
%       up = up0;
%       return;
%    end  
% %%
%    alp = 1; alpconst = 0.5; 
%    for iter = 1:maxit
%       xi = xi0 + alp*dxi;
%       Atxi = Atxi0+alp*Atdxi;
%       uinput = x0 - sig*Atxi;
%       [up,info_u] = proxFL(Binput,uinput,sig*lambda1,sig*lambda2);
%       par.info_u = info_u;
%       u = (uinput -up)/sig;
%       Ly = -b'*xi - 0.5*norm(xi)^2  - norm(up)^2/(2*sig);
%       if printlevel
%          fprintf('\n iter = %3d line search ------------------------------------- \n',iter);
%          fprintf('\n alp = %4.3f, LQ = %11.10e, LQ0 = %11.10e',alp,Ly,Ly0);
%          fprintf('\n ------------------------------------- \n');
%       end
%       if Ly-Ly0-c1*alp*g0 > -1e-12/max(1,abs(Ly0))
%          if (printlevel); fprintf(':'); end
%          return
%       else
%          alp = alpconst*alp;
% %          if printlevel; 
% %             fprintf('\n iter = %2d, ------line search value------------\n',iter);
% %             fprintf('\n ------alp = %2.2f, LQ = %11.10e, LQ0 = %11.10e',alp,Ly,Ly0);       
% %          end
%       end
% %      keyboard
%    end 
% %   if iter == maxit; keyboard; end
%    if (printlevel); fprintf('m'); end
%********************************************************************
function [par,Ly,xi,Atxi,u,up,alp,iter] = ...
         findstep(par,Binput,b,lambda1,lambda2,Ly0,xi0,Atxi0,...
         u0,up0,dxi,Atdxi,tol,options)

   if isfield(options,'stepop'); stepop = options.stepop; end
   printlevel =0; 
   maxit = ceil(log(1/(tol+eps))/log(2));
   c1 = 1e-4; c2 = 0.9; 
   sig = par.sigma;
  % Amap = par.lsAmap;
%%  
   g0  = dxi'*(-b - xi0) + Atdxi'*up0;
   %Ly0 = -b'*xi0 - 0.5*norm(xi0)^2  - norm(up0)^2/(2*sig);
   %grad0 = -(xi0 + b - Amap(up0));
   if (g0 <= 0)
      alp = 0; iter = 0; 
      if (printlevel)
         fprintf('\n Need an ascent direction, %2.1e  ',g0); 
      end
      xi = xi0;
      Atxi = Atxi0;
      u = u0;
      up = up0;
      Ly = Ly0;
      return;
   end  
%%
   alp = 1; alpconst = 0.5; 
   proxop = 2;
   for iter = 1:maxit
      if (iter==1)         
         alp = 1; LB = 0; UB = 1; 
      else
         alp = alpconst*(LB+UB);
      end
      xi = xi0 + alp*dxi;
      %Atxi = Atxi0+alp*Atdxi;
      %uinput = x0 - sig*Atxi;
      uinput = up0 + sig*u0 - sig*alp*Atdxi;
      if proxop == 1
         [up,info_u] = proxFL(Binput,uinput,sig*lambda1,sig*lambda2);
         par.info_u = info_u;
      else
         utmp = mexCondat(uinput,sig*lambda2);
         up = sign(utmp).*max(abs(utmp) - sig*lambda1,0);
      end
      u = (uinput -up)/sig;
      galp = dxi'*(-b - xi) + Atdxi'*up;
      Ly = -b'*xi - 0.5*norm(xi)^2  - norm(up)^2/(2*sig);
      %gradalp = -(xi + b - Amap(up));
      if printlevel
          fprintf('\n ------------------------------------- \n');
          fprintf('\n alp = %4.3f, LQ = %11.10e, LQ0 = %11.10e',alp,Ly,Ly0);
          fprintf('\n galp = %4.3f, g0 = %4.3f',galp,g0);
          fprintf('\n ------------------------------------- \n');
      end
      if (iter==1)
         gLB = g0; gUB = galp; 
         if (sign(gLB)*sign(gUB) > 0)
            if (printlevel); fprintf('|'); end
            Atxi = Atxi0+alp*Atdxi;
            if proxop ~=1
               info_u.rr1 = (abs(up) > 0);
               info_u.rr2 = (abs(Binput.Bmap(utmp)) < 1e-12);
               par.info_u = info_u;
            end
            return;             
         end
      end
      if ((abs(galp) < c2*abs(g0)))         
         if (Ly-Ly0-c1*alp*g0 > -1e-12/max(1,abs(Ly0))) ...
            && (stepop==1) || ((stepop == 2) && (abs(galp) < tol))
            if (printlevel); fprintf(':'); end
            Atxi = Atxi0+alp*Atdxi;
            if proxop ~=1
               info_u.rr1 = (abs(up) > 0);
               info_u.rr2 = (abs(Binput.Bmap(utmp)) < 1e-12);
               par.info_u = info_u;
            end
            return
         end
      end
      if (sign(galp)*sign(gUB) < 0)
         LB = alp; gLB = galp;
      elseif (sign(galp)*sign(gLB) < 0) 
         UB = alp; gUB = galp; 
      end
   end 
   if (iter == maxit)
      Atxi = Atxi0+alp*Atdxi;
      if proxop ~=1
         info_u.rr1 = (abs(up) > 0);
         info_u.rr2 = (abs(Binput.Bmap(utmp)) < 1e-12);
         par.info_u = info_u;
      end
   end
   if (printlevel); fprintf('m'); end
% %%********************************************************************
         

