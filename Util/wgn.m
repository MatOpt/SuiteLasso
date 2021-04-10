function y = wgn(varargin)
%WGN Generate white Gaussian noise.
%   Y = WGN(M,N,P) generates an M-by-N matrix of white Gaussian noise. P
%   specifies the power of the output noise in dBW. The unit of measure for
%   the output of the wgn function is Volts. For power calculations, it is
%   assumed that there is a load of 1 Ohm. 
%
%   Y = WGN(M,N,P,IMP) specifies the load impedance in Ohms.
%
%   Y = WGN(M,N,P,IMP,S) uses S, which is a random stream handle, to
%   generate random noise samples with RANDN. This syntax is useful to
%   generate repeatable outputs.  Type 'help RandStream' for more
%   information.
%
%   Y = WGN(M,N,P,IMP,STATE) resets the state of RANDN to STATE. This usage
%   is deprecated and may be removed in a future release.
%
%   Additional flags that can follow the numeric arguments are:
%
%   Y = WGN(..., POWERTYPE) specifies the units of P.  POWERTYPE can be
%   'dBW', 'dBm' or 'linear'.  Linear power is in Watts.
%
%   Y = WGN(..., OUTPUTTYPE); Specifies the output type.  OUTPUTTYPE can be
%   'real' or 'complex'.  If the output type is complex, then P is divided
%   equally between the real and imaginary components.
%
%   Example 1: 
%       % To generate a 1024-by-1 vector of complex noise with power
%       % of 5 dBm across a 50 Ohm load, use:
%       Y = wgn(1024, 1, 5, 50, 'dBm', 'complex')
%
%   Example 2: 
%       % To generate a 256-by-5 matrix of real noise with power
%       % of 10 dBW across a 1 Ohm load, use:
%       Y = wgn(256, 5, 10, 'real')
%
%   Example 3: 
%       % To generate a 1-by-10 vector of complex noise with power
%       % of 3 Watts across a 75 Ohm load, use:
%       Y = wgn(1, 10, 3, 75, 'linear', 'complex')
%
%   See also RANDN, AWGN.

%   Copyright 1996-2014 The MathWorks, Inc.

% --- Initial checks
narginchk(3,7);

% --- Value set indicators (used for the strings)
pModeSet    = 0;
cplxModeSet = 0;

% --- Set default values
pMode    = 'dbw';
imp      = 1;
cplxMode = 'real';
seed     = [];

% --- Placeholders for the numeric and string index values
numArg = [];
strArg = [];

% --- Identify string and numeric arguments
%     An empty in position 4 (Impedance) or 5 (Seed) are considered numeric
isStream = false;
for n=1:nargin
   if(isempty(varargin{n}))
      switch n
      case 4
         if(ischar(varargin{n}))
            error(message('comm:wgn:InvalidDefaultImp'));
         end;
         varargin{n} = imp; % Impedance has a default value
      case 5
         if(ischar(varargin{n}))
            error(message('comm:wgn:InvalidNumericInput'));
         end;
         varargin{n} = [];  % Seed has no default
      otherwise
         varargin{n} = '';
      end;
   end;

   % --- Assign the string and numeric vectors
   if(ischar(varargin{n}))
      strArg(size(strArg,2)+1) = n; %#ok<AGROW>
   elseif(isnumeric(varargin{n}))
      numArg(size(numArg,2)+1) = n; %#ok<AGROW>
   elseif(isa(varargin{n},'RandStream'))
      numArg(size(numArg,2)+1) = n; %#ok<AGROW>
      isStream = true;
   else
      error(message('comm:wgn:InvalidArg'));
   end;
end;

% --- Build the numeric argument set
switch(length(numArg))

   case 3
      % --- row is first (element 1), col (element 2), p (element 3)

      if(all(numArg == [1 2 3]))
         row    = varargin{numArg(1)};
         col    = varargin{numArg(2)};
         p      = varargin{numArg(3)};
      else
         error(message('comm:wgn:InvalidSyntax'))
      end;

   case 4
      % --- row is first (element 1), col (element 2), p (element 3), imp (element 4)
      %

      if(all(numArg(1:3) == [1 2 3]))
         row    = varargin{numArg(1)};
         col    = varargin{numArg(2)};
         p      = varargin{numArg(3)};
         imp    = varargin{numArg(4)};
      else
         error(message('comm:wgn:InvalidSyntax'))
      end;

   case 5
      % --- row is first (element 1), col (element 2), p (element 3), imp (element 4), seed (element 5)

      if(all(numArg(1:3) == [1 2 3]))
         row    = varargin{numArg(1)};
         col    = varargin{numArg(2)};
         p      = varargin{numArg(3)};
         imp    = varargin{numArg(4)};
         seed   = varargin{numArg(5)};
      else
         error(message('comm:wgn:InvalidSyntax'));
      end;
   otherwise
      error(message('comm:wgn:InvalidSyntax'));
end;

% --- Build the string argument set
for n=1:length(strArg)
   switch lower(varargin{strArg(n)})
   case {'dbw' 'dbm' 'linear'}
      if(~pModeSet)
         pModeSet = 1;
         pMode = lower(varargin{strArg(n)});
      else
         error(message('comm:wgn:TooManyPowerTypes'));
      end;
   case {'db'}
      error(message('comm:wgn:InvalidPowerType'));
   case {'real' 'complex'}
      if(~cplxModeSet)
         cplxModeSet = 1;
         cplxMode = lower(varargin{strArg(n)});
      else
         error(message('comm:wgn:TooManyOutputTypes'));
      end;
   otherwise
      error(message('comm:wgn:InvalidArgOption'));
   end;
end;

% --- Arguments and defaults have all been set, either to their defaults or by the values passed in
%     so, perform range and type checks

% --- p
if(isempty(p))
   error(message('comm:wgn:InvalidPowerVal'));
end;

if(any([~isreal(p) (length(p)>1) (isempty(p))]))
   error(message('comm:wgn:InvalidPowerVal'));
end;

if(strcmp(pMode,'linear'))
   if(p<0)
      error(message('comm:wgn:NegativePower'));
   end;
end;

% --- Dimensions
if(any([isempty(row) isempty(col) ~isscalar(row) ~isscalar(col)]))
   error(message('comm:wgn:InvalidDims'));
end;

if(any([(row<=0) (col<=0) ~isreal(row) ~isreal(col) ((row-floor(row))~=0) ((col-floor(col))~=0)]))
   error(message('comm:wgn:InvalidDims'));
end;

% --- Impedance
if(any([~isreal(imp) (length(imp)>1) (isempty(imp)) any(imp<=0)]))
   error(message('comm:wgn:InvalidImp'));
end;

% --- Seed
if(~isempty(seed))
    if ~isStream
        validateattributes(seed, {'double', 'RandStream'}, ...
            {'real', 'scalar', 'integer'}, 'WGN', 'S');
    end
end;

% --- All parameters are valid, so no extra checking is required
switch lower(pMode)
   case 'linear'
      noisePower = p;
   case 'dbw'
      noisePower = 10^(p/10);
   case 'dbm'
      noisePower = 10^((p-30)/10);
end;

% --- Generate the noise
if(~isempty(seed))
    if isStream
        hStream = seed;
    else
        hStream = RandStream('shr3cong', 'Seed', seed);
    end
    func = @(a,b)randn(hStream,a,b);
else
    func = @(a,b)randn(a,b);
end;

if(strcmp(cplxMode,'complex'))
   y = (sqrt(imp*noisePower/2))*(func(row, col)+1i*func(row, col));
else
   y = (sqrt(imp*noisePower))*func(row, col);
end;
