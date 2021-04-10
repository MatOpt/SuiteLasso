 function ys = ifftn_fast(xs)
%function ys = ifftn_fast(xs)
%
% For some reason, matlab's ifftn routine is suboptimal
% for the case of 2D FFTs, at least on some machines.
% The improvement herein was found by Hugo Shi.
%
% Note: matlab's ifft() and ifftn() handle an optional second "N" argument in
% different ways!  So to be safe I am not allowing any second argument here.
%
% Copyright 2004-6-28, Jeff Fessler, The University of Michigan

if nargin ~= 1, help(mfilename), error(mfilename), return, end
if streq(xs, 'test'), ifftn_fast_test, return, end

if ndims(xs) == 2 % 2D or 1D cases
	if min(size(xs)) == 1 % 1D
		ys = ifft(xs);
	else
		ys = ifftn_fast_fftfft(xs);
	end
else
	ys = ifftn(xs);
end


function ys = ifftn_fast_fftfft(xs)
ys = ifft(ifft(xs).').';


% test configuration of ifftn_fast for this machine
function ifftn_fast_test
ifftn_fast_test2
%ifftn_fast_test3


% test configuration of ifftn_fast for this machine for 2D
function ifftn_fast_test2

n = 2^8;
x = rand(n,n);

printf('starting test; be patient.')

% first loop is to get everything in cache or whatever.
% doing it twice is the only way to get an accurate comparison!
for nloop = [2 40];
	tic, for ii=1:nloop, y{1} = ifftn_fast(x); end
	tt(1) = toc; ty{1} = 'fftn_fast';

	tic, for ii=1:nloop, y{2} = ifftn(x); end
	tt(2) = toc; ty{2} = 'fftn';

	tic, for ii=1:nloop, y{3} = ifft(ifft(x).').'; end
	tt(3) = toc; ty{3} = 'fftfft_inline';

	tic, for ii=1:nloop, y{4} = ifft(ifft(x, [], 1), [], 2); end
	tt(4) = toc; ty{4} = 'fftfft_brack';

	tic, for ii=1:nloop, y{5} = ifft2(x); end
	tt(5) = toc; ty{5} = 'ifft2';

	tic, for ii=1:nloop, y{6} = ifftn_fast_fftfft(x); end
	tt(6) = toc; ty{6} = 'fftfft_func';
end

for ii = 1:length(tt)
	printf('time %14s = %g', ty{ii}, tt(ii))
	if max_percent_diff(y{1}, y{ii}) > 1e-11, error 'bug', end
end

if tt(1) > 1.10 * min(tt(2:end))
	error 'ifftn_fast is configured supoptimally for your machine!'
else
	printf('ifftn_fast is configured appropriately for your machine')
	printf('ifftn / ifftn_fast = %g%% ', tt(2) / tt(1) * 100.)
end
