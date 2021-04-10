function v = mypartition(n, L1)
% Chose (n-1) the splitting points of the array [0:(n+L1)]
s = nchoosek(1:n+L1-1,n-1);
m = size(s,1);

s1 = zeros(m,1,class(L1));
s2 = (n+L1)+s1;


v = diff([s1 s s2],1,2); 
v = v-1;

end 