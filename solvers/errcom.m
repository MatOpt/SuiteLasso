  function etaorg = errcom(x,grad,lambda1,lambda2)
      tmp = x - grad;
      tmp1 = mexCondat(tmp,lambda2);
      tmp2 = sign(tmp1).*max(abs(tmp1) - lambda1,0);
      etaorg = norm(x - tmp2);