function L = mychol(M,m)
%% 1/13
%% chol decompostion for M
   pertdiag = 1e-15*ones(m,1); 
   M = M + spdiags(pertdiag,0,m,m);
%%
   if (nnz(M) < 0.2*m*m); use_spchol=1; else; use_spchol=0; end
   if (use_spchol)
      [L.R,L.p,L.perm] = chol(sparse(M),'vector'); 
      L.Rt = L.R'; 
      L.matfct_options = 'spcholmatlab'; 
   else
      if issparse(M); M = full(M); end;           
      L.matfct_options = 'chol';
      L.perm = [1:m]; 
      [L.R,indef] = chol(M); 
   end
%%***********************************************************