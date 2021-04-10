%%******************************************************
%% Run this script in Matlab command window to generate some mex files
%% 
%% LASSO: 
%% Copyright (c) 2017 by
%% Xudong Li, Defeng Sun, and Kim-Chuan Toh 
%%******************************************************

   function Installmex(recompile) 

   if (nargin == 0); recompile = 0; end
   computer_model = computer;
   matlabversion = sscanf(version,'%f');
   matlabversion = matlabversion(1);
   tmp = version('-release'); 
   matlabrelease = str2num(tmp(1:4));
%%
   mexcmd = 'mex -O  -largeArrayDims  -output ';    
%%
   if (matlabversion < 7.3)  && (matlabrelease <= 2008)
      error(' needs MATLAB version 7.4 and above'); 
   end
   fsp = filesep;
   libstr = [];   

   curdir = pwd;  
   fprintf(' current directory is:  %s\n',curdir); 
%%
%% generate mex files in mexfun
%%
   clear fname

   src = [curdir,fsp,'mexfun']; 
   eval(['cd ','mexfun']); 
   fprintf ('\n Now compiling the mexFunctions in:\n'); 
   fprintf (' %s\n',src);       
   %%
   fname{1} = 'mexbwsolve';
   fname{2} = 'mexfwsolve';
   fname{3} = 'mextriang'; 
   fname{4} = 'mexMatvec';   
   fname{5} = 'mexscale';
   fname{6} = 'mexsigma_update_Classic_Lasso_SSNAL';   
   fname{7} = 'mexsigma_update_Fused_Lasso_SSNAL';
   fname{8} = 'mexCondat';
   fname{9} = 'mexFusedLassoJacobian';

   existfile = zeros(1,length(fname));
   for k = 1:length(fname)
      existfile(k) = exist([fname{k},'.',mexext]);   
   end
   for k = 1:length(fname)
       if recompile || ~existfile(k)
          cmd([mexcmd,fname{k},'  ',fname{k},'.c  ',libstr]);  
       end
   end      
   fprintf ('\n Compilation of mexFunctions completed.\n'); 
   cd .. 
   %cd ..
%%***********************************************
   function cmd(s) 
   
   fprintf(' %s\n',s); 
   eval(s); 
%%***********************************************
