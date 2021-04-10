%%**********************************************************************
%%**********************************************************************
    function [xi,resnrm,solve_ok] = Classic_Lasso_linsys_solver(Ainput,rhs,par)
    
    m = length(rhs);
    pp = ~par.rr;
    Ayes = isfield(Ainput,'A');
    solver = 'd_pcg';
    dn = 10000;
    sp = sum(pp);
%%
    if (m <= dn) && Ayes  
       if (m <= 1000) %1000
          solver = 'd_direct';
       elseif sp <= max(0.01*par.n,dn)
          solver = 'd_direct';
       end
    end
    if sp <= 0.7*m && Ayes && sp <=dn 
       solver = 'p_direct';
    end
    if (m > 5e3 && sp >= 200) || (m>2000 && sp > 800) || (m > 100 && sp > 1e4)
       solver = 'd_pcg';
    end
%%
    if strcmp(solver,'d_pcg')
       if Ayes
          AP = Ainput.A(:,pp);
          if false
             tmp = sum(AP.*AP,2);
             par.precond = 1;
             par.invdiagM = 1./(1 + par.sigma*tmp);
          end
          [xi,~,resnrm,solve_ok] = ...
          psqmry('matvec_ClassicLasso',AP,rhs,par); 
       else
          [xi,~,resnrm,solve_ok] = ...
          psqmry('matvec_ClassicLasso_Amap',Ainput,rhs,par); 
       end
    elseif strcmp(solver,'d_direct')
       AP = Ainput.A(:,pp);
       sigAPAt = par.sigma*(AP*AP');
       if m <= 1500      
          M = eye(m) + sigAPAt;
          xi = M\rhs;
       else
          M = speye(m,m) + sigAPAt;  
          L = mychol(M,m);
          xi = mylinsysolve(L,rhs);
       end
       resnrm = 0; solve_ok = 1;
    elseif strcmp(solver,'p_direct')
       AP = Ainput.A(:,pp);
       APT = AP';
       rhstmp = APT*rhs; 
       PAtAP = APT*AP;
       if sp <= 1500
          M = eye(sp)/par.sigma + PAtAP;
          tmp = M\rhstmp;
       else
          M = speye(sp,sp)/par.sigma + PAtAP;
          L = mychol(M,sp);
          tmp = mylinsysolve(L,rhstmp);
       end
       resnrm = 0; solve_ok = 1;
       xi = rhs - AP*tmp;
    end
%%**********************************************************************
