 %%*********************************************************************
 %%*********************************************************************
   clear all
   rng('default');
   addpath(genpath(pwd));
   warning off
   %%
   HOME = pwd; 
   addpath(genpath(HOME));
    
  fname{2}  = 'blocksig';
  fname{9}  = 'blkheavi';
  fname{10} = 'blknheavi';
  fname{701} = 'blurrycam';
  fname{702} = 'blurspike';
  fname{3} = 'cosspike';
  fname{11} = 'gausspike';
  fname{5}  = 'gcosspike';
  fname{902} = 'jitter';
  fname{6} = 'p3poly';
  fname{7} = 'sgnspike';
  fname{903} = 'spiketrn';
  fname{601} = 'soccer1';
  fname{602} = 'soccer2';
  fname{401} = 'srcsep1';
  fname{402} = 'srcsep2';
  fname{403} = 'srcsep3';
  fname{603} = 'yinyang';
  
  fid  = [2,3,5,6,7,9:11,401:403,601:603, 701, 702, 902,903];
  noiseyes = 1;
%%
for i = [403] 
    fprintf('\n Problem name: %s', fname{i});
    prob = generateProblem(i);
    Amap = @(x) prob.A(x,1);
    ATmap = @(x) prob.A(x,2);
    a = prob.sizeA;
    m = a(1); n = a(2);
    b = prob.b;
    if noiseyes
       dR = 60;
       b = awgn(b,dR,'measured');
       ns = 'ns';
    else
       ns = [];
    end
    AATmap = @(x) Amap(ATmap(x));
    eigsopt.issym = 1;
    Lip = eigs(AATmap,length(b),1,'LA',eigsopt);
    fprintf('\n Lip const = %3.2e, normb= %3.2e \n', Lip, norm(b));
    fprintf('\n-----------------------------------------------');
    fprintf('------------------------------')
    fprintf('\n     norm(ATy,inf) = %g',norm(ATmap(b),inf));

    for crho = [3] 
       rho = (10^(-crho))*max(abs(ATmap(b)));
       fprintf('\n Parameters: rho = %g',rho);
       fprintf('\n-----------------------------------------------');
       fprintf('------------------------------')
       fprintf('\n Problem ID = %3.0d, Problem size: n = %g,  m = %g',i,n,m)
       fprintf('\n-----------------------------------------------');
       fprintf('------------------------------')
       stoptol = 1e-6;
       if (true)
          opts.stoptol = stoptol;
          opts.Lip = Lip;
          Ainput.Amap = @(x) Amap(x);
          Ainput.ATmap = @(x) ATmap(x);
          [obj,y,xi,x,info,runhist] = Classic_Lasso_SSNAL(Ainput,b,n,rho,opts);
          Snal_res = info;
       end
    end
end
%%*********************************************************************

