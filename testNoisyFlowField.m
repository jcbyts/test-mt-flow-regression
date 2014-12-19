%% random vector field
% Hyperflow is too correlated over space and time. We can recover spatial
% RFs, but we have no resolution to resolve time. Important reminder: If
% you have a lot of correlated data, you effectively have very little data.
%

% stimulus dimensions and duration
if ~exist('ny', 'var') % y dimension of sitmulus
    ny = 8;
end

if ~exist('nx', 'var') % x dimension of sitmulus
    nx = 8;
end

if ~exist('nt', 'var') % number of frames in stimulus
    nt = 60e3;
end

if ~exist('ntk', 'var') % number of time points to smooth with
    ntk = 20;
end

if ~exist('temporalSmoothness', 'var') % time constant of smoothness (frames)
    temporalSmoothness   = 2;
end

if ~exist('spatialSmoothness', 'var') % spatial smoothness
    spatialSmoothness = 3;
end

if ~exist('spatres', 'var')
    spatres = 1;
end

% make noise
nse = randn(ny*nx*2+ny+1, nt);

% get stimulus dimension
xax = linspace(-nx/2,nx/2,nx);
yax = linspace(-ny/2,ny/2,ny);

[xx,yy] = meshgrid(xax,yax);

% Build filter to smooth noise.. motion requires some correlation in time
% and space to be motion, so we'll lowpass white noise using a
% spatiotemporal kernel
spatF = exp(-(xx.^2 + yy.^2)/spatialSmoothness); % spatial filter is a 2D gaussian

tr    = 0:ntk; % time range

spatK = spatF(:); % spatial kernel
timF  = exp(-(tr-ntk/2).^2/temporalSmoothness); % temporal filter is also gaussian
%timF  = tr.*exp(-tr/tau)/tau; % temporal filter is an alpha function

stK = spatK*timF; % spatiotemporal kernel

% Figure 1: noise and smoothing kernel
figure(1); clf

subplot(131) % space
imagesc(xax, yax, spatF); axis xy
xlabel('x space')
xlabel('y space')
title('spatial smoothing kernel')

subplot(132) % time
plot(tr, timF, 'k')
xlabel('time (frames)')
ylabel('amplitude')
title('temporal kernel')

subplot(133) % space-time
imagesc(stK)
xlabel('time')
ylabel('2D space')
title('full spatiotemporal smoothing kernel')
drawnow
disp('Create the stimulus by generting noise in velocity space')
disp('Smoothing white noise will generate correlations at the spatial and temporal scale of the smoothing kernel')
commandwindow

% convolve smoothing kernel with noise
nfilt = conv2(nse, stK, 'same');

% figure 2, filtered noise
figure(2); clf
subplot(231)
imagesc(nse(:,1:500)');
title('white noise sample')
subplot(232)
imagesc(fftshift(real(fft2(nse'))))
title('fft of noise')
subplot(233)
X = nfilt(ceil(nx/2)+1:end-ceil(nx/2)-1,:)';
imagesc(cov(nse'))
title('covariance matrix of white noise')
subplot(234)
imagesc(X(1:500,:));
title('filtered noise sample')
subplot(235)
imagesc(fftshift(real(fft2(X))))
title('fft of filtered noise')
subplot(236)
imagesc(cov(X))
title('covariance of filtered noise')
drawnow
disp('A look at the white noise once it has been filtered shows that there')
disp('are plenty of spatiotemporal correlations')

clear nfilt nse