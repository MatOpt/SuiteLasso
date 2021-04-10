%%*****************************************************************
    function [up,info] = proxFL(Binput,u, lambda1, lambda2)
    op = 1;
    utmp = mexCondat(u,lambda2);
    if (nargout > 1)
       [up,rr1] = proxL1(utmp,lambda1);
       if op == 1
          tmp = Binput.Bmap(utmp); 
          rr2 = (abs(tmp) < 1e-12);
       else
          ztmp = mexRose(Binput.Bmap(u - utmp));
          [~,rr2] = proj_inf(ztmp,lambda2);
       end
       info.rr1 = rr1;
       info.rr2 = rr2;
       info.innerNT = 0;
       info.innerflsa = 0;
    else
       up = sign(utmp).*max(abs(utmp) - lambda1,0);
    end
%%*****************************************************************

