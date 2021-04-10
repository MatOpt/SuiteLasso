
function [y,Atxi,xi,par,runhist,info] = ...
    Classic_Lasso_SSNCG(n,b,Ainput,x0,Ax0,Atxi0,xi0,ld,par,options)

printsub = 1;
breakyes = 0;
maxitersub = 50;
tiny = 1e-10;
tol = 1e-6;
maxitpsqmr =500;
precond = 0;
Ascaleyes = 0;

if isfield(options,'printsub'); printsub = options.printsub; end
if isfield(options,'maxitersub'); maxitersub = options.maxitersub; end
if isfield(options,'tiny'); tiny = options.tiny; end
if isfield(options,'tol'); tol = options.tol; end
if isfield(options,'precond'); precond = options.precond; end
if isfield(options,'Ascaleyes'); Ascaleyes = options.Ascaleyes; end

existA = options.existA;
sig = par.sigma;
bscale = options.bscale;
cscale = options.cscale;
normborg = 1+norm(b)*sqrt(bscale*cscale);
%% preperation
Amap = @(x) Ainput.Amap(x);
ATmap = @(x) Ainput.ATmap(x);
yinput = -Atxi0 - x0/sig;
[y,rr] = proj_inf(yinput,ld);
par.rr = rr;
Rpb = Ax0 - b + xi0;
normRp = norm(Rpb);
ytmp = yinput - y;
Atxi = Atxi0; xi = xi0;
Ly = b'*xi - 0.5*norm(xi)^2 - 0.5*sig*norm(ytmp)^2;
runhist.psqmr(1) = 0;
runhist.findstep(1) = 0;
cnt_Amap = 0;
cnt_ATmap = 0;
cnt_pAATmap = 0;
cnt_fAATmap = 0;
%% mian Newton iteration
for itersub = 1:maxitersub    
    yold = y; xiold = xi; 
    Atxiold = Atxi; 
    Rdz = Atxi + y;
    normRd = norm(Rdz);
    %Ly = b'*xi - 0.5*norm(xi)^2 - 0.5*sig*norm(ytmp)^2;
    msigAytmp = -sig*Amap(ytmp);
    GradLxi = -(xi - b + msigAytmp); %-(xi - b -sig*Amap(ytmp));
    cnt_Amap = cnt_Amap + 1;
    normGradLxi = norm(GradLxi)*sqrt(bscale*cscale)/normborg;
    priminf_sub = normGradLxi; 
    if Ascaleyes == 1
       dualinf_sub = norm(Rdz./options.dscale)*cscale/(1+norm(y./options.dscale)*cscale);
    else
       dualinf_sub = normRd*cscale/(1+norm(y)*cscale);
    end
    if max(priminf_sub,dualinf_sub) < tol
       tolsubconst = 0.9;
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
    if max([normGradLxi]) < max(tolsub) && itersub > 1
        msg = 'good termination in subproblem:';
        if printsub
            fprintf('\n       %s  ',msg);
            fprintf(' dualinfes = %3.2e, gradLyxi = %3.2e, tolsub = %3.2e',...
                          dualinf_sub,normGradLxi,tolsub);
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
     par.epsilon = min(1e-3,0.1*normGradLxi); %% good to add
     par.precond = precond;
     if precond == 1
        par.invdiagM = 1/(1+sig);
     end
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
     if Ascaleyes == 1 && false
        tolpsqmr = min(5e-3, 0.1*norm(rhs));
     else
        tolpsqmr = min(5e-3, 0.1*norm(rhs));
     end
     const2 = 1;
     if itersub > 1 && (prim_ratio > 0.5 || priminf_sub > 0.1*runhist.priminf(1))
        const2 = 0.5*const2;
     end
     if (dual_ratio > 1.1); const2 = 0.5*const2; end
     tolpsqmr = const2*tolpsqmr;
     par.tol = tolpsqmr; par.maxit = maxitpsqmr;
     [dxi,resnrm,solve_ok] = Classic_Lasso_linsys_solver(Ainput,rhs,par);
     Atdxi = ATmap(dxi);
     cnt_ATmap = cnt_ATmap + 1;
     iterpsqmr = length(resnrm)-1;
     if (iterpsqmr ==0)
        cnt_pAATmap = cnt_pAATmap + 1; 
     else
         if existA
            cnt_pAATmap = cnt_pAATmap + iterpsqmr;
         else
            cnt_fAATmap = cnt_fAATmap + iterpsqmr; 
         end
     end
     if (printsub)
          fprintf('| %3.1e %3.1e %3.0d %-3d',par.tol,resnrm(end),iterpsqmr);
          fprintf(' %2.1f %2.0d',const2,sum(1-par.rr));
     end
     par.iter = itersub;
     if (itersub<=3) && (dualinf_sub > 1e-4) || (par.iter <3)
         stepop = 1;
     else
         stepop = 2;
     end
     steptol = 1e-5; step_op.stepop=stepop;
     [par,Ly,xi,Atxi,y,ytmp,alp,iterstep] = ...
         findstep(par,b,ld,Ly,xi,Atxi,y,ytmp,dxi,Atdxi,steptol,step_op); 
     runhist.solve_ok(itersub) = solve_ok;
     runhist.psqmr(itersub)    = iterpsqmr; 
     runhist.findstep(itersub) = iterstep; 
     if alp < tiny; breakyes =11; end
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
             Rdz =  Atxi + y;
             msigAytmp = -sig*Amap(ytmp);
             cnt_Amap = cnt_Amap + 1;
             %normRd = norm(Rdz);
             if printsub
                if Ascaleyes
                   fprintf('\n new dualfeasorg = %3.2e', norm(Rdz./options.dscale)*cscale/(1+norm(y./options.dscale)*cscale));
                else
                   fprintf('\n new dualfeasorg = %3.2e', norm(Rdz)*cscale/(1+norm(y)*cscale));
                end
             end
             break
          end
     end 
end
info.maxCG = max(runhist.psqmr);
info.avgCG = sum(runhist.psqmr)/itersub;
info.breakyes = breakyes;
info.itersub = itersub;
info.tolconst = par.tolconst;
info.RpGradratio  = normRp*sqrt(bscale*cscale)/(normGradLxi*normborg);
info.rankX = par.rr;
info.ytmp = ytmp;
info.cnt_Amap = cnt_Amap;
info.cnt_ATmap = cnt_ATmap;
info.Ax = msigAytmp;
info.cnt_pAATmap = cnt_pAATmap;
info.cnt_fAATmap = cnt_fAATmap;
%%********************************************************************
%%********************************************************************
function [par,Ly,xi,Atxi,y,ytmp,alp,iter] = ...
         findstep(par,b,ld,Ly0,xi0,Atxi0,y0,ytmp0,dxi,Atdxi,tol,options)
     
   if isfield(options,'stepop'); stepop = options.stepop; end
   printlevel = 0; 
   maxit = ceil(log(1/(tol+eps))/log(2));
   c1 = 1e-4; c2 = 0.9; 
   sig = par.sigma;
%%
   tmp1 = dxi'*(b-xi0);
   tmp2 = norm(dxi)^2; 
   g0  = tmp1 + sig*Atdxi'*ytmp0; 
   Ly = []; 
   if (g0 <= 0)
      alp = 0; iter = 0; 
      if (printlevel) 
         fprintf('\n Need an ascent direction, %2.1e  ',g0); 
      end
      xi = xi0;
      Atxi = Atxi0;
      y = y0;
      ytmp = ytmp0;
      Ly = Ly0;
      return;
   end  
%%
   alp = 1; alpconst = 0.5; 
   for iter = 1:maxit
      if (iter==1)          
         alp = 1; LB = 0; UB = 1; 
      else
         alp = alpconst*(LB+UB);
      end
      xi = xi0 + alp*dxi;
      yinput = ytmp0 + y0 - alp*Atdxi;
      [y,rr] = proj_inf(yinput,ld);
      par.rr = rr;
      ytmp = yinput - y;
      galp = tmp1 - alp*tmp2 + sig*Atdxi'*ytmp;
      if (iter==1)
         gLB = g0; gUB = galp; 
         if (sign(gLB)*sign(gUB) > 0)
            if (printlevel); fprintf('|'); end
            Atxi = Atxi0+alp*Atdxi;
            Ly = b'*xi - 0.5*norm(xi)^2 - 0.5*sig*norm(ytmp)^2;             
            return;             
         end
      end
      if (abs(galp) < c2*abs(g0)) 
         Ly = b'*xi - 0.5*norm(xi)^2 - 0.5*sig*norm(ytmp)^2;
         if (Ly-Ly0-c1*alp*g0 > -1e-8/max(1,abs(Ly0))) ...
            && ((stepop==1) || (stepop==2 && abs(galp)<tol))
            if (printlevel); fprintf(':'); end
            Atxi = Atxi0+alp*Atdxi;
            return;            
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
   end
   if (printlevel); fprintf('m'); end
   if isempty(Ly)
      Ly = b'*xi - 0.5*norm(xi)^2 - 0.5*sig*norm(ytmp)^2;            
   end
%%********************************************************************
         

