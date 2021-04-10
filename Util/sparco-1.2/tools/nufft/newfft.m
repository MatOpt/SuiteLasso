 function st = newfft(om_in, Nd_in, varargin)
%function st = nufft_init(om, Nd, [options])
%
% UNDER CONSTRUCTION
% New version of NUFFT (pun intended) that uses real interpolation kernels.
% (The original NUFFT code used complex interpolation needlessly.)
%
% This returns a "strum" object with methods for both forward and adjoint
% d-dimensional NUFFT operations.  \sum_{n=0}^{N-1} x[n] exp(-1i * om * n)
%
% in
%	om [M,d]	"digital" frequencies in radians (can be empty!)
%				(if empty, then user provides 'om' to methods)
%	Nd [d]		image dimensions (N1,N2,...,Nd)
%
% options
%	'Jd' [d]	# of neighbors used (in each direction). (default: 6)
%	'Kd' [d]	FFT sizes (should be >= Nd). (default: 2*Nd)
%	n_shift [d]	n = 0-n_shift to N-1-n_shift (default: 0)
%		Like fft(), the NUFFT expects the signals to be x(0,0), ...
%		Use n_shift = [N1/2, N2/2, ...] for x(-N1/2,-N2/2,...), ...
%
%	'mode'	char	how to compute the NUFFT (default: 'table0')
%			'table0' table with nearest neighbor interpolation
%			'table1' table with linear interpolation
%			'sparse' precompute large sparse matrix
%				caution: this option may require lots of memory!
%
%	'ktype'	char	what type of interpolation kernel (default: 'minmax:kb')
%
% ktype options:
%	'exact'		slow FT
%	'kb:beatty'	KB with parameters from Beatty et al T-MI Jun 2005
%	'minmax:kb'	minmax interpolator with excellent KB scaling!
%	'minmax:tuned'	minmax interpolator, somewhat numerically tuned
%	'minmax:user'	minmax interpolator with user ({alpha}, {beta})
%	'uniform'	uniform scaling factors (not recommended)
%	'kaiser'	kaiser-bessel (KB) interpolator (minmax best alpha, m)
%			or 'kaiser', [alpha], [m] to specify parameter (vectors)
%	'linear'	linear interpolator (a terrible straw man)
%	kernel		user-provided inline interpolation kernel(k,J)
%			(or a cell array of kernels, one for each dimension)
%	'table'		use table-based interpolation rather than sparse matrix.
%			this can save a lot of memory for large problems.
%			example ..., 'table', 2^11, 'minmax:kb'
%			where 2^11 is the table over-sampling factor.
% out
%	st	strum object with several methods including the following:
%		st.fft(x, [om])
%		st.adj(Xo, [om])
%		st.p		[M, *Kd]	sparse interpolation matrix
%						(or empty if table-based) 
%		st.sn		[(Nd)]		scaling factors
%		st.Nd,Jd,Kd,om	copies of inputs
%
% *Nd is shorthand for prod(Nd).
% (Nd) is shorthand for (N1,N2,...,Nd)
%
% Copyright 2007-6-3, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(om_in, 'test'), newfft_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% inputs
st.om = om_in; % [M,d] frequency samples
st.Nd = Nd_in; % [d] signal dimentions
st.dd = length(st.Nd); % dimensionality of input space (usually 2 or 3)

% defaults
st.Jd = 6 * ones(1, st.dd);
st.Kd = 2 * st.Nd;
st.n_shift = zeros(1, st.dd);
st.mode = 'table0';
st.ktype = 'kaiser';
st.alpha = {}; % minmax parameters
st.beta = {};
st.tol = 0;
st.is_kaiser_scale = false;

% options
st = vararg_pair(st, varargin);

% special cases of input sampling pattern
if ischar(st.om)
	st.om = nufft_samples(st.om, st.Nd);
end

% checks
if st.dd ~= length(st.Jd) | st.dd ~= length(st.Kd)
	error 'inconsistent dim'
end
if st.dd ~= length(st.n_shift)
	error('n_shift needs %d columns', st.dd)
end

if any(st.Kd < st.Nd), warning 'Kd < Nd unlikely to work.  Try Kd=2*Nd', end

if ~isempty(st.om) && st.dd ~= size(st.om,2)
	error('omega needs %d columns', st.dd)
end

if streq(st.mode, 'sparse') && isempty(st.om)
	error 'sparse mode requires "om"'
end

%
% special cases
%
if streq(st.ktype, 'exact')
	st = strum(st, { ...
		'fft', @newfft_exact_fft, '(x, [om])';
		'adj', @newfft_exact_adj, '(X, [om])';
		'p', @newfft_exact_p, '([om])';
		'sn', @(s) ones(s.Nd), '()';
		});
return
end

error 'not done'

%
% process optional arguments
%

%switch st.ktype
%case

% table based?
if ischar(varargin{1}) & streq(varargin{1}, 'table')
	st = nufft_table_init(om, Nd, Jd, Kd, n_shift, varargin{2:end});
	return
end

ktype = varargin{1};

% cell array of kernel functions: {kernel1, kernel2, ..., kernelD}
if isa(ktype, 'cell')
	if isa(ktype{1}, 'inline') | isa(ktype{1}, 'function_handle')
		ktype = 'inline';
		if length(varargin) > 1, error 'excess arguments?', end
		if length(varargin{1}) ~= dd, error 'wrong # of kernels', end
		st.kernel = varargin{1};
	else
		error 'cell array should be inline kernels!?'
	end

% or a single inline kernel for all dimension
elseif isa(ktype, 'inline') | isa(ktype, 'function_handle')
	ktype = 'inline';
	if length(varargin) > 1, error 'excess arguments?', end
	for id = 1:dd
		st.kernel{id} = varargin{1};	% all same
	end

% or a string that describes the type of interpolator
elseif ~ischar(ktype)
	error 'non-string kernel type?'

end

%
% set up whatever is needed for each interpolator
%
if streq(ktype, 'inline')
	% already did it above

% linear interpolator straw man
elseif streq(ktype, 'linear')
	ktype = 'inline';
	kernel = inline('(1 - abs(k/(J/2))) .* (abs(k) < J/2)', 'k', 'J');
	for id = 1:dd
		st.kernel{id} = kernel;
	end

% KB interpolator
elseif streq(ktype, 'kaiser')
	is_kaiser_scale = true;

	% with minmax-optimized parameters
	if length(varargin) == 1
		for id = 1:dd
			[st.kernel{id} st.kb_alf(id) st.kb_m(id)] = ...
				kaiser_bessel('inline', Jd(id));
		end

	% with user-defined parameters
	elseif length(varargin) == 3
		alpha_list = varargin{2};
		m_list = varargin{3};
		if (length(alpha_list) ~= dd) | (length(m_list) ~= dd)
			error('#alpha=%d #m=%d vs dd=%d', ...
				length(alpha_list), length(m_list), dd)
		end
		for id = 1:dd
			[st.kernel{id} st.kb_alf(id) st.kb_m(id)] = ...
				kaiser_bessel('inline', Jd(id), ...
					alpha_list(id), m_list(id));
		end
	else
		error 'kaiser should have no arguments, or both alpha and m'
	end

% minmax interpolator with KB scaling factors (recommended default)
elseif streq(ktype, 'minmax:kb')
	for id = 1:dd
		[st.alpha{id}, st.beta{id}] = ...
			nufft_alpha_kb_fit(Nd(id), Jd(id), Kd(id));
	end

% minmax interpolator with numerically "tuned" scaling factors
elseif streq(ktype, 'minmax:tuned')
	for id = 1:dd
		[st.alpha{id}, st.beta{id}, ok] = ...
			nufft_best_alpha(Jd(id), 0, Kd(id)/Nd(id));
		if ~ok, error 'unknown J,K/N', end
	end

% minmax interpolator with user-provided scaling factors
elseif streq(ktype, 'minmax:user')
	if length(varargin) ~= 3, error 'user must provide alpha/beta', end
	st.alpha = varargin{2};
	st.beta = varargin{3};
	if length(st.alpha) ~= dd | length(st.beta) ~= dd
		error 'alpha/beta size mismatch'
	end

elseif streq(ktype, 'uniform')
	for id = 1:dd
		st.alpha{id} = 1;
		st.beta{id} = 0;
	end

else
	error 'unknown kernel type'
end

st.M = size(om,1);

%
% scaling factors: "outer product" of 1D vectors
%
st.sn = 1;
for id=1:dd
	if is_kaiser_scale
		nc = [0:Nd(id)-1]'-(Nd(id)-1)/2;
		tmp = 1 ./ kaiser_bessel_ft(...
			nc/Kd(id), Jd(id), st.kb_alf(id), st.kb_m(id), 1);
	elseif streq(ktype, 'inline')
		tmp = 1 ./ nufft_interp_zn(0, Nd(id), Jd(id), Kd(id), st.kernel{id});
	else
		tmp = nufft_scale(Nd(id), Kd(id), st.alpha{id}, st.beta{id});
	end
	st.sn = st.sn(:) * tmp';
end
if length(Nd) > 1
	st.sn = reshape(st.sn, Nd);	% [(Nd)]
else
	st.sn = st.sn(:);	% [(Nd)]
end

%
% [J?,M] interpolation coefficient vectors.  will need kron of these later
%
for id=1:dd
	N = Nd(id);
	J = Jd(id);
	K = Kd(id);
	if isvar('st.kernel')
		[c, arg] = ...
		nufft_coef(om(:,id), J, K, st.kernel{id});	% [J?,M]
	else
		alpha = st.alpha{id};
		beta = st.beta{id};
		T = nufft_T(N, J, K, st.tol, alpha, beta);	% [J?,J?]
		[r, arg] = ...
		nufft_r(om(:,id), N, J, K, alpha, beta);	% [J?,M]
		c = T * r;	clear T r
	end

	gam = 2*pi/K;
	phase_scale = 1i * gam * (N-1)/2;

	phase = exp(phase_scale * arg);	% [J?,M] linear phase
	ud{id} = phase .* c;		% [J?,M]

	%
	% indices into oversampled FFT components
	%
	koff = nufft_offset(om(:,id), J, K);	% [M,1] to leftmost near nbr
	kd{id} = mod(outer_sum([1:J]', koff'), K) + 1;	% [J?,M] {1,...,K?}
	if id > 1	% trick: pre-convert these indices into offsets!
		kd{id} = (kd{id}-1) * prod(Kd(1:(id-1)));
	end

end, clear c arg gam phase phase_scale koff N J K

%
% build sparse matrix that is [M,*Kd]
% with *Jd nonzero entries per frequency point
%
if dd >= 3, printm('Needs at least %g Gbyte RAM', prod(Jd)*M*8/10^9*2), end

kk = kd{1};	% [J1,M]
uu = ud{1};	% [J1,M]
for id = 2:dd
	Jprod = prod(Jd(1:id));
	kk = block_outer_sum(kk, kd{id});	% outer sum of indices
	kk = reshape(kk, Jprod, M);
	uu = block_outer_prod(uu, ud{id});	% outer product of coefficients
	uu = reshape(uu, Jprod, M);
end	% now kk and uu are [*Jd, M]

%
% apply phase shift
% pre-do Hermitian transpose of interpolation coefficients
%
phase = exp(1i * (om * n_shift(:))).';			% [1,M]
uu = conj(uu) .* phase(ones(1,prod(Jd)),:);		% [*Jd,M]

mm = [1:M]; mm = mm(ones(prod(Jd),1),:);		% [*Jd,M]
st.p = sparse(mm(:), kk(:), uu(:), M, prod(Kd));	% sparse matrix
% sparse object, to better handle single precision operations!
st.p = Gsparse(st.p, 'odim', [M 1], 'idim', [prod(Kd) 1]);


%
% in
%	x1	[J1,M]
%	x2	[J2,M]
% out
%	y	[J1,J2,M]	y(i1,i2,m) = x1(i1,m) + x2(i2,m)
%
function y = block_outer_sum(x1, x2)
[J1 M] = size(x1);
[J2 M] = size(x2);
xx1 = reshape(x1, [J1 1 M]);	% [J1,1,M] from [J1,M]
xx1 = xx1(:,ones(J2,1),:);	% [J1,J2,M], emulating ndgrid
xx2 = reshape(x2, [1 J2 M]);	% [1,J2,M] from [J2,M]
xx2 = xx2(ones(J1,1),:,:);	% [J1,J2,M], emulating ndgrid
y = xx1 + xx2;			% [J1,J2,M]

function y = block_outer_prod(x1, x2)
[J1 M] = size(x1);
[J2 M] = size(x2);
xx1 = reshape(x1, [J1 1 M]);	% [J1,1,M] from [J1,M]
xx1 = xx1(:,ones(J2,1),:);	% [J1,J2,M], emulating ndgrid
xx2 = reshape(x2, [1 J2 M]);	% [1,J2,M] from [J2,M]
xx2 = xx2(ones(J1,1),:,:);	% [J1,J2,M], emulating ndgrid
y = xx1 .* xx2;			% [J1,J2,M]


%
% newfft_test()
%
function newfft_test
Nd = [20 10];
st_e = newfft('epi', Nd, 'ktype', 'exact');
rand('state', 0)
x0 = rand(Nd);

% test exact
Xe = st_e.fft(x0);
xe = st_e.adj(Xe);
pe = st_e.p;
equivs(pe * x0(:), Xe)

'todo'
