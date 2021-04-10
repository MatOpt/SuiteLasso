%%*********************************************************************
%%*********************************************************************
clear all
rng('default');
%%
HOME = pwd; 
addpath(genpath(HOME));
datadir = [HOME,filesep,'UCIdata']; 

fname{1} = 'E2006.train';
fname{2} = 'log1p.E2006.train';
fname{3} = 'E2006.test';
fname{4} = 'log1p.E2006.test';
fname{5} = 'pyrim_scale_expanded5';
fname{6} = 'triazines_scale_expanded4';
fname{7} = 'abalone_scale_expanded7';
fname{8} = 'bodyfat_scale_expanded7';
fname{9} = 'housing_scale_expanded7';
fname{10} = 'mpg_scale_expanded7';
fname{11} = 'space_ga_scale_expanded9';

%%
for i = [5]
   probname = [datadir,filesep,fname{i}];
   fprintf('\n Problem name: %s', fname{i});
   if exist([probname,'.mat'])
      load([probname,'.mat'])
   else
      fprintf('\n Can not find the file in UCIdata');
      fprintf('\n ');
      return
   end   
   [m,n] = size(A);
   if exist('mexMatvec')
      Amap = @(x) mexMatvec(A,x,0);
      ATmap = @(y) mexMatvec(A,y,1);
   else
      AT = A';
      Amap  = @(x) A*x;
      ATmap = @(x) AT*x; 
   end
   AATmap = @(x) Amap(ATmap(x));
   eigsopt.issym = 1;
   Lip = eigs(AATmap,length(b),1,'LA',eigsopt);
   fprintf('\n Lip const = %3.2e, nomrb = %3.2e ', Lip, norm(b));
   for crho = 4;[3,4]; 
       c = 10^(-crho);
       rho = c*max(abs(ATmap(b)));
       fprintf('\n--------------------------------------------------');
       fprintf('------------------------------')
       fprintf('\n       rho = %g *norm(ATy,inf)',c);
       fprintf('\n Parameters: rho = %g',rho);
       fprintf('\n--------------------------------------------------');
       fprintf('------------------------------')
       fprintf('\n Problem: n = %g,  m = %g',n,m)
       fprintf('\n--------------------------------------------------');
       fprintf('------------------------------')
       %% 
       stoptol = 1e-6;
       if (true)
          opts.stoptol = stoptol;
          opts.Lip = Lip;
          opts.Ascale = 1;
          Ainput.A = A;
          Ainput.Amap = @(x) Amap(x);
          Ainput.ATmap = @(x) ATmap(x);
          [obj,y,xi,x,info,runhist] = Classic_Lasso_SSNAL(Ainput,b,n,rho,opts);
          Snal_res = info;
       end
   end
end
%%*********************************************************************

