clear all
rng('default');

filepath = mfilename('fullpath');
idxpwd = find(filepath == filesep);
lidxpwd = length(idxpwd);

foldpath = filepath(1:idxpwd(lidxpwd)-1);
addpath(genpath(foldpath));

srcdatadir = [foldpath,filesep,'UCIdataorg'];

HOME = filepath(1:idxpwd(lidxpwd-2)-1);
newdatadir = [HOME,filesep,'UCIdata']; 

if ~exist(newdatadir,'dir')
   mkdir(newdatadir);
end

expandconst = ones(11,1);
fname{1} = 'E2006.train';
fname{2} = 'log1p.E2006.train';
fname{3} = 'E2006.test';
fname{4} = 'log1p.E2006.test';
fname{5} = 'pyrim_scale';     expandconst(5) = 5;
fname{6} = 'triazines_scale'; expandconst(6) = 4;
fname{7} = 'abalone_scale';   expandconst(7) = 7;
fname{8} = 'bodyfat_scale';   expandconst(8) = 7;
fname{9} = 'housing_scale';   expandconst(9) = 7;
fname{10} = 'mpg_scale';      expandconst(10) = 7;
fname{11} = 'space_ga_scale'; expandconst(11) = 9;


for i = 8; %1:11
   probname = fname{i};
   if i>= 5; probnamenew = [probname,'.txt']; else probnamenew = probname; end
   if exist(probnamenew)
      [b,A] = libsvmread([srcdatadir,filesep,probnamenew]);
   else
      fprintf('\n i = %d, prob %s, no such file in UCIdataorg!', i, probname);
      return
   end
   [m,n] = size(A);
   %%
   % expand dataset
   if i>= 5
      if i == 5; fprintf('\n '); end
      fprintf('\n problem %s, m = %3.0d, n = %3.0d ', fname{i},...
          m,n);
      d = expandconst(i);
      v = mypartition(n+1,d);
      v = v';
      nn = size(v,2);
      fprintf('\n after expansion, n ---> %3.0d', nn);
      AA = zeros(m,nn);
      for q=1:nn
         AAq = ones(m,1);
         for j = 1:n
            if v(j,q) > 0
              AAq = AAq.*(A(:,j).^v(j,q));
            end
         end
         AA(:,q) = AAq;
      end
      A = AA;
      probname = [probname,'_expanded',num2str(d)];
      n = nn;
   end
   
   %%
   % remove zero columns
   aa = full(sqrt(sum(A.*A)));
   idx = find(aa > 0);
   if length(idx) < n
      A = A(:,idx);
      fprintf('\n Problem %s removed %d zero columns', fname{i}, n - length(idx));
   end
   
   %%
   % save data
   save([newdatadir,filesep,probname,'.mat'],'b','A');
   fprintf('\n Done! problem ID = %3.0d, problemname = %s ', i,...
       fname{i});
   clear b A probname n
end
fprintf('\n ');
 