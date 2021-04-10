    function q = mylinsysolve(L,r) 

%     if (L.isidentity)
%        q = r; 
%     else
       if strcmp(L.matfct_options,'chol')
          q(L.perm,1) = mextriang(L.R, mextriang(L.R,r(L.perm),2) ,1);
       elseif strcmp(L.matfct_options,'spcholmatlab')
          q(L.perm,1) = mexbwsolve(L.Rt,mexfwsolve(L.R,r(L.perm,1)));
       end
    end