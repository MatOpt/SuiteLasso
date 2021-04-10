%%****************************************************************
%% solve I + sigma*A*J*At
%%****************************************************************

function invM = FL_precond_directinverse(Ainput,par)


   rr1 = par.info_u.rr1; %%for L1
   rr2 = par.info_u.rr2; %%for fused  
   [h,U] = Jacobian(rr2,rr1); 
   if ~isempty(U)
      PU = U(rr1,:);
   else
      PU = []; PU1 = 0;
   end
   Ph = h(rr1);    
   nrm = sum(PU); 
   nzcolidx = find(nrm > 0); 
   numblk1 = length(nzcolidx);    
   if (numblk1==0) && (norm(Ph)==0)
      invM = @(rhs) rhs; 
      return
   end
   if (numblk1 > 0)
      PU1 = PU(:,nzcolidx);          
   end
   AP = Ainput.A(:,rr1);  

 %% direct solver                  
   m  = size(Ainput.A,1);
   lenP = length(Ph);        
   V1 = AP(:,Ph==1); 
   if (nnz(V1) > 0.1*numel(V1)) && issparse(V1)
      V1 = full(V1); 
   end
   if (numblk1 < 2*lenP)
      if (m < lenP)
          V2 = AP*PU1; 
          if (nnz(V2) > 0.1*numel(V2)) && issparse(V2)
             V2 = full(V2); 
          end             
          M = speye(m,m) + par.sigma*(V1*V1'+V2*V2');     
          L = mychol(M,length(M));      
          invM = @(rhs) mylinsysolve(L,rhs);  
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
             invM = @(rhs) mylinsysolve(L,rhs);         
          else
             nW = size(W,2) ;   
             if (nnz(W) > 0.1*numel(W)) && issparse(W)
                W = full(W); 
             end
             SMWmat = W'*W;             
             SMWmat = spdiags(ones(nW,1)/par.sigma,0,nW,nW)+ SMWmat;     
             L = mychol(SMWmat,nW);
             invM = @(rhs) rhs - W*mylinsysolve(L,(rhs'*W)');       
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
       if (numblk1 > 0.5*m)          
          M = speye(m,m) + par.sigma*(V1*V1'+V2*AP');     
          L = mychol(M,length(M));      
          invM = @(rhs) mylinsysolve(L,rhs);         
       else
          nW = size(V1,2)+size(V2,2); 
          V = [V1,V2]; W = [V1,AP];
          SMWmat = spdiags(ones(nW,1)/par.sigma,0,nW,nW)+W'*V;
          if (nnz(SMWmat)/nW^2 < 0.4)
             [L.L,L.U,L.p,L.q] = lu(sparse(SMWmat),'vector');    
          else
             [L.L,L.U,L.p] = lu(full(SMWmat),'vector');
             L.q = [1:nW]; 
          end
          invM = @(rhs) invMfun(L,W,nW,rhs,V);
       end       
   end
end
%%**************************************************************
%% Jacobian = spdiags(hh,0,n,n) + U*U';
%%**************************************************************
   function [hh,U] = Jacobian(rr2,rr1)
        
    n = length(rr2)+1; 
    if (false)   
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
   end
%%****************************************************************
%% final invM
%%****************************************************************
   function xi = invMfun(L,W,nW,rhs,V)
      tmp = (rhs'*W)'; 
      tmp2 = zeros(nW,1);
      tmp2(L.q) = L.U\(L.L\tmp(L.p)); 
      xi = rhs - V*tmp2;         
   end
%%****************************************************************
