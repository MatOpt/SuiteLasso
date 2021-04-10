%%****************************************************************
%% solve I + sigma*A*J*At
%%****************************************************************

function [xi,resnrm,solve_ok,par] = Fused_Lasso_linsys_solver(Ainput,rhs,par)

   solver_options = 2; %%2=direct solver
   resnrm = 0; solve_ok = 1;
   rr1 = par.info_u.rr1; %%for L1
   rr2 = par.info_u.rr2; %%for fused  
   m = length(rhs);
   if (sum(rr1) > 5e3 && length(rhs) > 5e3) || (sum(rr1)>1.5e4 && length(rhs)>2000) || ~isfield(Ainput,'A') %% for large P use CG
      solver_options = 1;
      par.precond = 0;
   end
   if (solver_options==1)
      par.cg = [];
   end
   if false
      [h,U] = Jacobian(rr2,rr1); 
   else
      %%NEW: much faster
      [h,U] = mexFusedLassoJacobian(double(rr2)); 
   end
   if ~isempty(U)
      Ut = U'; 
      PU = Ut(:,rr1)';
   else
      PU = []; PU1 = 0;
   end
   Ph = h(rr1);    
   nrm = sum(PU); 
   nzcolidx = find(nrm > 0); 
   numblk1 = length(nzcolidx);
   lenP = length(Ph);
   par.lenP = lenP;
   par.numblk1 = numblk1;
   if (numblk1==0) && (norm(Ph)==0)
      xi = rhs; 
      par.innerop = 2;
      return; 
   end
   if (numblk1 > 0)
      PU1 = PU(:,nzcolidx);
   else
      PU1 = 0;
   end
   if isfield(Ainput,'A') 
      AP = Ainput.A(:,rr1);        
      if (solver_options==1)
      %% iterative solver
         par.AP = AP; 
         par.PU1 = PU1;
         par.Ph  = Ph;        
         if ~isfield(par.cg,'solve') par.cg.solve = 0; end
         if ~isfield(par.cg,'iter') par.cg.iter = 0; end
         if par.cg.solve ~= 1 || par.cg.iter > 50
            if m <= 6000 && sum(rr1) <= 1e4 %% using direct solver to compute preconditioner
               par.precond = 3;               
               par.invM = FLprecond_direct_inverse(Ainput,par);
            else
               [dA,VA,VAt] = FLprecond_partial_eigendecomp(m,20,par);
               par.dA = dA;
               par.V = VA;
               par.Vt = VAt;
               par.precond = 2;
            end
         end
        if isfield(par,'dA'); par.d = 1 + par.sigma*par.dA; end
         [xi,~,resnrm,solve_ok] = psqmry('matvec_FusedLasso',Ainput,rhs,par);
         par.cg.solve = solve_ok;
         par.cg.iter = length(resnrm);
      elseif (solver_options==2)
      %% direct solver                  
         m  = size(Ainput.A,1); 
         V1 = AP(:,Ph==1); 
         if (nnz(V1) > 0.1*numel(V1)) && issparse(V1)
                    V1 = full(V1); 
         end
         if (numblk1 < 2*lenP)
            if (m < lenP + numblk1)
               V2 = AP*PU1;        
               if (nnz(V2) > 0.1*numel(V2)) && issparse(V2)
                  V2 = full(V2); 
               end         
               tmpM1 = par.sigma*(V1*V1'+V2*V2');
               M = speye(m,m) + tmpM1;     
               L = mychol(M,length(M));      
               xi = mylinsysolve(L,rhs);         
            else
               if (numblk1 > 0)
                  W = [V1, AP*PU1];             
               else
                  W = V1;  
               end
               if (numblk1 > 0.5*m)      
                  if (nnz(W) > 0.1*numel(W)) && issparse(W)
                     W = full(W); 
                  end  
                  M = speye(m,m) + par.sigma*(W*W');     
                  L = mychol(M,length(M));      
                  xi = mylinsysolve(L,rhs);     
               else
                  nW = size(W,2) ;
                  if issparse(W)
                     if (nnz(W) > 0.1*numel(W)) %&& issparse(W)
                        W = full(W); 
                     end
                  end
                  SMWmat = W'*W;             
                  SMWmat = spdiags(ones(nW,1)/par.sigma,0,nW,nW)+ SMWmat; 
                  L = mychol(SMWmat,nW);
                  xi = rhs - W*mylinsysolve(L,(rhs'*W)');       
               end
            end
         else
            if (numblk1 > 0)
               V2 = AP*(PU1*PU1');
               if (nnz(V2) > 0.1*numel(V2)) && issparse(V2)
                  V2 = full(V2); 
               end 
            else 
               V2 = [];
            end        
            nW = size(V1,2)+size(V2,2); 
            if m < nW % (numblk1 > 0.5*m) %        
               M = speye(m,m) + par.sigma*(V1*V1'+V2*AP');     
               L = mychol(M,length(M));      
               xi = mylinsysolve(L,rhs);         
            else
               V = [V1,V2]; W = [V1,AP];
               SMWmat = spdiags(ones(nW,1)/par.sigma,0,nW,nW)+W'*V;
               if (nnz(SMWmat)/nW^2 < 0.4)
                  [L.L,L.U,L.p,L.q] = lu(sparse(SMWmat),'vector');    
               else
                  [L.L,L.U,L.p] = lu(full(SMWmat),'vector');
                  L.q = [1:nW]; 
               end
               tmp = (rhs'*W)'; 
               tmp2 = zeros(nW,1);
               tmp2(L.q) = L.U\(L.L\tmp(L.p));
               xi = rhs - V*tmp2;      
            end       
         end
      end
   else
      par.PU1 = PU1;
      par.Ph  = Ph;    
      [xi,~,resnrm,solve_ok] = psqmry('matvec_FusedLasso_Amap',Ainput,rhs,par);
      par.cg.solve = solve_ok;
      par.cg.iter = length(resnrm);
   end
  par.innerop = solver_options;
%%**************************************************************
%% Jacobian = spdiags(hh,0,n,n) + U*U';
%%**************************************************************
   function [hh,U] = Jacobian(rr2,rr1)
        
    n = length(rr2)+1; 
    if (true)   
       blklen =[]; blkend = [];
       len = 0; numblk = 0;            
       for k=1:length(rr2)
          if (rr2(k)==1) 
             len=len+1; 
          else
             if (len > 0)
                numblk = numblk+1;
                blklen(numblk,1) = len;
                blkend(numblk,1) = k; 
                len = 0;                    
             end
          end
       end
       if (len > 0)
          numblk = numblk+1;
          blklen(numblk,1) = len;
          blkend(numblk,1) = n; 
       end
    else 
       %% much faster
       idx1 = find(rr2==1);
       idxblk = find(diff(idx1) > 1);    
       blkend = [];
       blklen = diff([0; idxblk]); 
       if ~isempty(idxblk)
          blkend = idx1(idxblk)+1;
       end
       if (rr2(end)==1)
          blkend(end+1) = length(rr2)+1;
          if ~isempty(idxblk)
             blklen(end+1) = length(rr2)-idx1(idxblk(end))-1;
          else
             blklen(end+1) = length(rr2); 
          end
       end
    end   
%%    
    numblk = length(blklen);
    NZ = sum(blklen)+numblk+1;
    ii = zeros(NZ,1); jj = zeros(NZ,1); vv = zeros(NZ,1);   
    hh = zeros(n,1); 
    cnt = 0;       
    for k=1:numblk
       len = blklen(k)+1; invsqrtlen = 1/sqrt(len);  
       idxend = blkend(k); 
       idxsub = [idxend-blklen(k): idxend];
       hh(idxsub) = 1; 
       ii(cnt+[1:len]) = idxsub;      
       jj(cnt+[1:len]) = k; %%k*ones(len,1);
       vv(cnt+[1:len]) = invsqrtlen; %%ones(len,1)/sqrt(len);
       cnt = cnt + len;
    end
    ii(cnt+1) = n; jj(cnt+1) = numblk; vv(cnt+1) = 0;
    if (numblk > 0)
       U = spconvert([ii,jj,vv]);
    else
       U = []; 
    end
    hh = ones(n,1)-hh; 
%%************************************************************
    function APAxi = APAt(xi,par)
    AP = par.AP;   
    PU1 = par.PU1;
    if (size(PU1,2) > 0)
       tmp = (xi'*AP)';  
       tmp2 = (tmp'*PU1)';
       APAxi = AP*(par.Ph.*tmp + PU1*tmp2);
    elseif (norm(par.Ph) > 0)
       tmp = (xi'*AP)';      
       APAxi = AP*(par.Ph.*tmp);
    else
       APAxi = zeros(length(xi)); 
    end
%%**************************************************************
%% partial eigen-decomposition preconditioner 
%%**************************************************************
   function [dA,VA,VAt] = FLprecond_partial_eigendecomp(m,rA,par)
   eigsopt.issym = 1;
   tmprA = 20;
   APAtmap = @(xi) APAt(xi,par);
   [VA,dA,flagA] = eigs(APAtmap,m,tmprA,'LA',eigsopt);
   dA = diag(dA); rA = min(rA,sum(dA>0));
   proxA = min(20,rA);
   if (proxA < tmprA)
      dA = dA(1:proxA);
      VA = VA(:,1:proxA);
   end
   VAt = VA';
%%**************************************************************
%% partial eigen-decomposition of A*AT
%%**************************************************************
   function [dA,VA,VAt] = AAT_partial_eigendecomp(Ainput,m,rA)
   eigsopt.issym = 1;
   tmprA = 10;
   AATmap = @(x) Ainput.Amap(Ainput.ATmap(x));
   [VA,dA,flagA] = eigs(AATmap,m,tmprA,'LA',eigsopt);
   dA = diag(dA); rA = min(rA,sum(dA>0));
   proxA = min(10,rA);
   if proxA < tmprA
      dA = dA(1:proxA);
      VA = VA(:,1:proxA);
   end
   VAt = VA';
%%**************************************************************   