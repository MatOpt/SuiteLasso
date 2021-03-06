function data = prob401(varargin)
%PROB401  Source separation example 1.
%
%   PROB401 creates a problem structure.  The generated signal will
%   consist of two mixed audio signals each two seconds long.
%
%   The following optional arguments are supported:
%
%   PROB401('mixing',MIXING,flags) is the same as above, but with the
%   2 by 2 mixing matrix set to MIXING. The flags include the default
%   'show' and 'update' flags for plotting the problem illustrations.
%
%   Examples:
%   P = prob401;  % Creates the default 401 problem.
%
%   References:
%
%   [freesound] Database of Creative Commons licensed sounds,
%     http://freesound.iua.upf.edu/ 
%
%   See also GENERATEPROBLEM.
%
%MATLAB SPARCO Toolbox.

%   Copyright 2008, Rayan Saab, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: prob401.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{},{'mixing'});
mixing      = getOption(parm, 'mixing', [0.1375, 0.7250; ...
                                         0.9905, 0.6888]);
info.name   = 'srcsep1';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Set up the data
s1 = audioread(sprintf('%sprob401_Guitar.wav', opts.datapath));%wavread(sprintf('%sprob401_Guitar.wav', opts.datapath));%
s1 = s1(1:6:end); % Downsample from 48.0kHz to 8.0kHz
s2 = audioread(sprintf('%sprob401_Piano.wav', opts.datapath));%wavread(sprintf('%sprob401_Piano.wav', opts.datapath));
s2 = s2(1:5:end); % Downsample from 44.1kHz to 8.8kHz
% Extract part
s1 = s1(:); s2 = s2(:);
s2 = s2(5000:5000+length(s1)-1);
m  = length(s1);

% Mix the signals
mixture = [s1 s2] * mixing';

% Set up the problem
data.signal          = [s1, s2];
data.op.MixingMatrix = opMatrix(kron(mixing,speye(m,m)), ...
                                'Mixing matrix: kron(mixing,speye)');
data.op.DCT          = opDCT(512);
data.op.WindowedDCT  = opWindowedOp(m, data.op.DCT, ones(512,1), 256);
data.op.BlockDiag    = opBlockDiag([1,1],data.op.WindowedDCT);
data.B               = data.op.BlockDiag;
data.M               = data.op.MixingMatrix;
data.b               = mixture(:);
data                 = completeOps(data);

% Additional information
info.title           = 'Source separation example 1';
info.thumb           = 'figProblem401';
info.citations       = {'freesound'};
info.fig{1}.title    = 'Audio signal 1';
info.fig{1}.filename = 'figProblem401Audio1';
info.fig{2}.title    = 'Audio signal 2';
info.fig{2}.filename = 'figProblem401Audio2';
info.fig{3}.title    = 'Mixed audio signal 1';
info.fig{3}.filename = 'figProblem401Mixed1';
info.fig{4}.title    = 'Mixed audio signal 2';
info.fig{4}.filename = 'figProblem401Mixed2';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(data.signal(:,1),'b-'); xlim([1,16000]);
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(data.signal(:,2),'b-'); xlim([1,16000]);
  updateFigure(opts, info.fig{2}.title, info.fig{2}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(data.b(1:m),'b-'); ylim([-1,1]); xlim([1,16000]);
  updateFigure(opts, info.fig{3}.title, info.fig{3}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(data.b(m+(1:m)),'b-'); ylim([-1,1]); xlim([1,16000]);
  updateFigure(opts, info.fig{4}.title, info.fig{4}.filename)

  if opts.update
     P1 = ones(128,128,3);
     P1 = thumbPlot(P1,1:m,data.signal(:,1),[0,0,1]);
     P2 = ones(128,128,3);
     P2 = thumbPlot(P2,1:m,data.signal(:,2),[1,0,0]);
     P  = 0.5 * P1 + 0.5 * P2;
     P = (P(1:2:end,:,:) + P(2:2:end,:,:)) / 2;
     P = (P(:,1:2:end,:) + P(:,2:2:end,:)) / 2;
     thumbwrite(P, info.thumb, opts);
  end
end
