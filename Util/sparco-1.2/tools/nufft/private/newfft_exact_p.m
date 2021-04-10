 function pre = newfft_exact_p(st, om)
% forward
nthread = 1; % todo
useloop = false; % todo

if ~isempty(st.om)
	om = st.om;
end

dd = st.dd;
Nd = st.Nd;
n_shift = st.n_shift;

if dd > 3
	error 'only up to 3D is done'
else
	Nd = [Nd(:); ones(3-length(Nd),1)];
	n_shift = [n_shift(:); zeros(3-length(n_shift),1)];
end

for id=1:3 % fix: dd
	nn{id} = [0:(Nd(id)-1)]-n_shift(id);
end

[nn{1} nn{2} nn{3}] = ndgrid(nn{1}, nn{2}, nn{3});

pre = 0;
for id=1:dd
	pre = pre + om(:,id) * col(nn{id})';
end
pre = exp(-1i*pre);
