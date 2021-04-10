%%*********************************************************************
%%*********************************************************************
clear all
rng('default');
%%
HOME = pwd; 
addpath(genpath(HOME));
%%
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

fname{21} = 'DLBCL_H';
fname{22} = 'lung_H1';
fname{23} = 'NervousSystem';
fname{24} = 'ovarian_P';
fname{25} = 'DLBCL_N';
fname{26} = 'DLBCL_S';
fname{27} = 'lung_H2';
fname{28} = 'lung_M';
fname{29} = 'lung_O';
fname{30} = 'ovarian_S';
%%
for i = [5] %[11];
   if i < 20
      datadir = [HOME,filesep,'UCIdata'];
   elseif i <= 30
      datadir = [HOME,filesep,'BioNUS'];
   end
   probname = [datadir,filesep,fname{i}];
   fprintf('\n Problem name: %s', fname{i});
   if exist([probname,'.mat'])
      load([probname,'.mat']) 
   else
      fprintf('\n SSNAL: can not find the file in UCIdata');
      fprintf('\n ');
      return
   end
   [m,n] = size(A);
   Amap  = @(x) A*x;
   ATmap = @(x) A'*x; 
   AATmap = @(x) Amap(ATmap(x));
   eigsopt.issym = 1;
   Lip = eigs(AATmap,length(b),1,'LA',eigsopt);
   fprintf('\n Lip const = %3.2e, nomrb = %3.2e ', Lip, norm(b));

   for crho = 3; %[7];
      c = 10^(-crho);
      lambda1 = c*max(abs(ATmap(b)));
      c2 = [10, 2, 0.2, 0.01, 1, 0.5];
   
      fprintf('\n-----------------------------------------------');
      fprintf('------------------------------')
      fprintf('\n Problem: n = %g,  m = %g',n,m)
      fprintf('\n-----------------------------------------------');
      fprintf('------------------------------')
   %%
      stoptol = 1e-6;

      for jj = 3; %[3,5,6];
         lambda2 = c2(jj)*lambda1;
         fprintf('\n lambda1 = %g *norm(ATy,inf)',c);
         fprintf('\n Parameters: lambda1 = %3.2e, lambda2 = %3.2e',lambda1, lambda2);
         opts.Lip = Lip;
         opts.stoptol = stoptol;
         opts.printyes = 1;
         inputformat = 1; 
         if (inputformat==1)
             Ainput = A; %%default
         elseif (inputformat==2)
             Ainput.A = A;
         elseif (inputformat==3)
             Ainput.Amap = @(x) A*x;
             Ainput.ATmap = @(x) A'*x;
         end
         [obj,x,xi,u,info,runhist] = ...
              Fused_Lasso_SSNAL(Ainput,b,n,lambda1,lambda2,opts);
      end
   end
end