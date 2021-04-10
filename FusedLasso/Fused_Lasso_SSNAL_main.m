function [obj,u,xi,x,info,runhist] = ...
              Fused_Lasso_SSNAL_main(Amap0,ATmap0,Bmap0,BTmap0,b,lambda1,lambda2,parmain,u,xi,x) 
   
   m = parmain.m; n = parmain.n;
   tstart  = parmain.tstart;
   scale   = parmain.scale;
   maxiter = parmain.maxiter;
   printyes= parmain.printyes;
   rescale = parmain.rescale;
   stoptol = parmain.stoptol;
   stopop  = parmain.stopop;
   precond = parmain.precond;
   Lip     = parmain.Lip;
   orgojbconst = parmain.orgojbconst;
   
   if isfield(parmain,'A'); A = parmain.A; end

   lambda1org = lambda1;
   lambda2org = lambda2;
   
   sigmaLip = 1/sqrt(Lip);
   Bmap = Bmap0;
   BTmap = BTmap0;
   borg = b;
   normborg = 1 + norm(borg);
   Atxi = ATmap0(xi); Ax = Amap0(x); Bx = Bmap(x);
   obj(1) = 0.5*norm(Ax - borg)^2 + lambda1org*norm(x,1) + lambda2org*norm(Bx,1) + orgojbconst;
   obj(2) = -(0.5*norm(xi)^2 + borg'*xi) + orgojbconst;
   bscale = 1; cscale = 1; 
   %%
   if (scale == 1)
      b = b/sqrt(bscale*cscale);
      xi = xi/sqrt(bscale*cscale);
      Amap = @(x) Amap0(x*sqrt(bscale/cscale));
      ATmap = @(x) ATmap0(x*sqrt(bscale/cscale));
      if exist('A','var'); A = A*sqrt(bscale/cscale); end
      lambda1 = lambda1/cscale; lambda2 = lambda2/cscale;
      x = x/bscale; u = u/cscale; 
      Bx = Bx/bscale;
      Ax = Ax/sqrt(bscale*cscale);
      Atxi = Atxi/cscale; 
      normb = 1+norm(b);
   end
   Ainput_nal.Amap = Amap;
   Ainput_nal.ATmap = ATmap;
   Binput_nal.Bmap = Bmap;
   Binput_nal.BTmap = BTmap;
   if exist('A','var'); Ainput_nal.A = A; end
   sigma = min(1, sigmaLip);
   if isfield(parmain,'sigma'); sigma = parmain.sigma; end
   Rp1 = Ax - b;
   Rd  = Atxi + u;
   normu = norm(u);
   normRp = norm(Rp1 - xi);
   normRd = norm(Rd);
   primfeas = normRp/normb;
   dualfeas = normRd/(1+normu);
   maxfeas = max(primfeas,dualfeas);
   dualfeasorg = normRd*cscale/(1+normu*cscale);
   primfeasorg = sqrt(bscale*cscale)*normRp/normborg;
   maxfeasorg = max(primfeasorg, dualfeasorg);
   runhist.dualfeasorg(1) = dualfeasorg;
   relgap = (obj(1) - obj(2))/(1+obj(1)+obj(2));
   if printyes
        fprintf('\n*******************************************************');
        fprintf('******************************');
        fprintf('\n \t\t   Fused Lasso: SSNAL ');
        fprintf('\n*******************************************************');
        fprintf('******************************');        
        fprintf('\n lambda1 = %3.2e, lambda2 = %3.2e',lambda1,lambda2);        
        fprintf('\n n = %3.0f, m = %3.0f',n, m);
        fprintf('\n bscale = %3.2e, cscale = %3.2e', bscale, cscale);
        fprintf('\n ---------------------------------------------------');
        fprintf('\n  iter|  [pinfeas  dinfeas]  [pinforg  dinforg]    relgaporg|    pobj          dobj    |');
        fprintf(' time | sigma |rankS |');
        fprintf('\n*****************************************************');
        fprintf('**************************************************************');
        fprintf('\n #%3.1d|  %3.2e %3.2e %3.2e %3.2e %- 3.2e %- 8.7e %- 8.7e  %5.1f',...
           0,primfeas,dualfeas,primfeasorg,dualfeasorg,relgap,obj(1),obj(2),etime(clock,tstart)); 
        fprintf('  %3.2e ',sigma);
   end
 %% ssncg
   SSNCG = 1;
   if SSNCG
      parNCG.sigma = sigma;
      parNCG.tolconst = 0.5;
      parNCG.n = n;
      parNCG.precond = precond;
   end
   maxitersub = 10;
   breakyes = 0;
   prim_win = 0;
   dual_win = 0;
   ssncgop.tol = stoptol;
   ssncgop.precond = precond;
   ssncgop.bscale = bscale;
   ssncgop.cscale = cscale;
   ssncgop.printsub = printyes;
   sigmamax = 1e8; sigmamin = 1e-4; 
   for iter = 1:maxiter
      if ((rescale == 1) && (maxfeas < 5e2) && (rem(iter,3) == 1) && (iter > 1) )...
         || ((rescale >= 2) && maxfeas < 1e-1 && (abs(relgap) < 0.05) ...
             && (iter >= 5) && (max(normx/normuxi,normuxi/normx) > 1.7) && rem(iter,5)==1)
         normAtxi = norm(Atxi);
         normx = norm(x); normu = norm(u);
         normuxi = max([normu,normAtxi]); 
         [sigma,bscale2,cscale2,sbc,sboc,bscale,cscale] = ...
                mexscale(sigma,normx,normuxi,bscale,cscale);
         b = b/sbc; 
         lambda1 = lambda1/cscale2;
         lambda2 = lambda2/cscale2;
         x = x/bscale2;
         xi = xi/sbc; 
         Atxi = Atxi/cscale2; Ax = Ax/sbc;
         u = u/cscale2; 
         Bx = Bx/bscale2;         
         Ainput_nal.Amap = @(x) Ainput_nal.Amap(x*sboc);
         Ainput_nal.ATmap = @(x) Ainput_nal.ATmap(x*sboc);
         if isfield(Ainput_nal,'A'); Ainput_nal.A = Ainput_nal.A*sboc; end
         if (precond == 2) && isfield(parNCG,'dA'); parNCG.dA = parNCG.dA*sboc; end
         ssncgop.bscale = bscale;
         ssncgop.cscale = cscale;
         normb = 1+norm(b);
         if printyes
            fprintf('\n    ');
            fprintf('[rescale=%1.0f: %2.0f| %3.2e %3.2e %3.2e | %3.2e %3.2e| %3.2e]',...
            rescale,iter,normx,normAtxi,normu,bscale,cscale,sigma);
         end         
         rescale = rescale+1; 
         prim_win = 0; dual_win = 0; 
      end 
      xold = x;  uold = u;
      Bxold = Bx;
      parNCG.sigma = sigma;
      parNCG.innerNT = 0;
      parNCG.innerflsa = 0;
      if dualfeas < 1e-5
         maxitersub = max(maxitersub,30);
      elseif dualfeas < 1e-3
         maxitersub = max(maxitersub,30);
      elseif dualfeas < 1e-1
         maxitersub = max(maxitersub,20);
      end
      ssncgop.maxitersub = maxitersub; 
      [u,Atxi,xi,parNCG,runhist_NCG,info_NCG] = ...
          Fused_Lasso_SSNCG(b,Ainput_nal,Binput_nal,x,Ax,Atxi,xi,...
                  lambda1,lambda2,parNCG,ssncgop);
      if info_NCG.breakyes < 0
         parNCG.tolconst = max(parNCG.tolconst/1.06,1e-3);
      end
      x = info_NCG.up;
      Ax = info_NCG.Aup;
      Rd = Atxi + u;
      %x = xold - gamma*sigma*Rd;
      %Ax = Ainput_nal.Amap(x);
      Bx = Bmap(x);
      Rp1 = Ax - b;
     %%----------------------------------------------------
      normRp = norm(Rp1 - xi);
      normRd = norm(Rd);
      normu = norm(u);
      primfeas = normRp/normb;
      dualfeas = normRd/(1+normu);
      maxfeas = max(primfeas,dualfeas);
      dualfeasorg = normRd*cscale/(1+normu*cscale);
      primfeasorg = sqrt(bscale*cscale)*normRp/normborg;
      maxfeasorg = max(primfeasorg, dualfeasorg);
      runhist.dualfeas(iter+1) = dualfeas;
      runhist.primfeas(iter+1) = primfeas;
      runhist.maxfeas(iter+1)  = maxfeas;
      runhist.primfeasorg(iter+1) = primfeasorg;
      runhist.dualfeasorg(iter+1) = dualfeasorg;
      runhist.maxfeasorg(iter+1)  = maxfeasorg;
      runhist.sigma(iter) = sigma;
      runhist.rank1(iter) = sum(parNCG.info_u.rr1);
      runhist.rank2(iter) = sum(parNCG.info_u.rr2);
      runhist.innerNT(iter) = parNCG.innerNT;
      runhist.innerflsa(iter) = parNCG.innerflsa;
      runhist.xr(iter) = sum(abs(x)>1e-10);
      runhist.Bxr(iter) = sum(abs(Bx)>1e-10);
%%---------------------------------------------------------
%% check for termination
%%--------------------------------------------------------
      if true
         objscale = cscale*bscale;
         primobj = objscale*(0.5*norm(xi)^2 + lambda1*norm(x,1) + lambda2*norm(Bx,1)) + orgojbconst;
         dualobj = objscale*(-0.5*norm(xi)^2 - b'*xi ) + orgojbconst;
         relgap = (primobj-dualobj)/( 1+abs(primobj)+abs(dualobj));
         ttime = etime(clock,tstart);
         %% check for termination
         if stopop == 1
            if (max([primfeasorg,dualfeasorg]) < 500*max(1e-6, stoptol)) 
                grad = ATmap0(Rp1*sqrt(bscale*cscale));
                etaorg = errcom(x*bscale,grad,lambda1org,lambda2org);
                %norm(x*bscale- proxFL(Binput,x*bscale - grad,lambda1org,lambda2org));
                eta = etaorg / (1 + norm(grad) + norm(x*bscale));
                if eta < stoptol 
                   breakyes = 1;
                   msg = 'converged';
                elseif abs(relgap) < stoptol && max([primfeasorg,dualfeasorg]) < stoptol && eta < sqrt(stoptol)
                   breakyes = 2;
                   msg = 'converged';
                end
            end
         elseif stopop == 2
            if max([primfeasorg,dualfeasorg]) < 2*stoptol
               grad = ATmap0(Rp1*sqrt(bscale*cscale));
               etaorg = errcom(x*bscale,grad,lambda1org,lambda2org);
               %etaorg = norm(x*bscale- proxFL(Binput,x*bscale - grad,lambda1org,lambda2org));
               eta = etaorg / (1 + norm(grad) + norm(x*bscale));
               if eta < stoptol
                  breakyes = 88;
                  msg = 'converged';
               end
            end
         end         
         if (printyes)
            fprintf('\n %5.0d| [%3.2e %3.2e] [%3.2e %3.2e]  %- 3.2e| %- 10.9e %- 10.9e |',...
               iter,primfeas,dualfeas,primfeasorg, dualfeasorg,relgap,primobj,dualobj); 
            fprintf(' %5.1f| %3.2e|',ttime, sigma); 
            if iter >= 1
               fprintf('%3d %3d|',sum(parNCG.info_u.rr1),sum(parNCG.info_u.rr2));
               fprintf('[%3d %3d]',sum(abs(x) > 1e-10), sum(abs(Bx) > 1e-10));
            end
            fprintf(' ttinnerNT = %3d ', sum(runhist.innerNT));
            if exist('eta'); fprintf('\n        [eta = %3.2e, etaorg = %3.2e]',eta, etaorg);end
         end
	     if (rem(iter,3*1)==1) 
            normx = norm(x); normAtxi = norm(Atxi); normu = norm(u);
            if (printyes)
               fprintf('\n        [normx,Atxi,u =%3.2e %3.2e %3.2e ]',...
               normx,normAtxi,normu);
            end
         end
         runhist.primobj(iter)   = primobj;
         runhist.dualobj(iter)   = dualobj;
         runhist.time(iter)      = ttime; 
         runhist.relgap(iter)    = relgap;
      end    
      if (breakyes > 0) 
         if printyes; fprintf('\n  breakyes = %3.1f, %s',breakyes,msg);  end
         break; 
      end
      if (primfeasorg < dualfeasorg)
         prim_win = prim_win+1; 
      else
         dual_win = dual_win+1; 
      end 
      [sigma,prim_win,dual_win] = ...
      mexsigma_update_Fused_Lasso_SSNAL(sigma,sigmamax,sigmamin,prim_win,dual_win,iter,parNCG.innerop);
   end
   
   if (iter == maxiter) && (breakyes == 0) 
      msg = 'maximum number of iterations reached';
   end
   
   info.iter = iter;
   info.bscale = bscale;
   info.cscale = cscale;
   info.Ax = Ax;
   info.ttime = ttime;
   info.msg = msg;
