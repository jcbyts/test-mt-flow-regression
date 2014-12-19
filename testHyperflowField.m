
if ~exist('frameRate', 'var')
    frameRate    = 60; % this is what we've got
end

if ~exist('nTotalFrames', 'var')
    nTotalFrames = 100e3;
end

if ~exist('ny', 'var') % y dimension of sitmulus
    ny = 8;
end

if ~exist('nx', 'var') % x dimension of sitmulus
    nx = 8;
end

if ~exist('spatres', 'var')
    spatres = 2;
end
   
if ~exist('maskdiameter', 'var')
    maskdiameter = nx/spatres;
end

if ~exist('filterOrder', 'var')
	filterOrder  = 4;
end

if ~exist('cutoffFreq', 'var')
    cutoffFreq   = 5;
end

% build params struct
params = struct();
params.designsizex = nx*spatres;
params.designsizey = ny*spatres;
params.maskdiameter = maskdiameter;
params.spatres = spatres; % spatial resolution

% build hyperflow inputs
fCutoff = cutoffFreq/(frameRate/2); % cutoff at 2 Hz
% we'll start with an order of 1 which should give us about 20 db/decade attenuation
[b,a] = butter(filterOrder,fCutoff); % use butter (because I know how... I think)

S = filter(b,a,randn(nTotalFrames, 6));

apertureCenter = randn(nTotalFrames+1,2)*params.designsizex; 
centerxy = filter(b,a,apertureCenter, [],1); 

X = GetVelField(params, S, centerxy);