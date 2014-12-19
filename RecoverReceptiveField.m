
%% Can we recover the neuron?
if ~exist('nPrincipalComponents', 'var')
    nPrincipalComponents = 1000;
end

% build design matrix
addDC = 1;
nk = size(X,2);
disp('building the design matrix')
B = makeStimRows(X, nkt, 'same'); % builds the design matrix
clear X
nY = numel(Y(:,1));
B = [B ones(nY, addDC)];
disp('precomputing useful quantities')
% pre compute some useful values
xx = B'*B; % stimulus covariance (unnormalized)
yy = Y'*Y; 
xy = B'*Y;
disp('running pca on the stimulus')
% stuff for principal components regression
[b,s] = svd(xx); % principal components of B
bas = b(:,1:nPrincipalComponents); % principal components basis

Xb = B*bas; % project the design matrix onto what the eigenvector basis
clear B % clear some memory for poor matlab

Xb = bsxfun(@minus, Xb, mean(Xb));

% do regression in the subspace
xxb = Xb'*Xb;
xyb = Xb'*Y;

wb  = xxb\xyb;
disp('done. Plot time')

%% plotting the Spike-Triggered Average
% the spike triggered average is easy to compute and it is proportional to
% the maximum-likelihood kernel if the stimulus is white noise. Ours isn't,
% so you'll see that the spatial and temporal RFs seem blured out by the
% correlations.

% this is the full spatiotemporal STA
STA = flipud(reshape(xy(1:end-addDC)/nY, nkt, nk));

% check if it's 
[u,~,v] = svd(STA);
n = 1;
sflip = sign(sum(u(1:floor(end/2),n)));

% plot the spike triggered average
fig = figure(14924); clf;
set(fig, 'Position', [1 515 1299 291]);
subplot(1, 3, 1)
x = v(1:end/2,n);
y = v(((end/2)+1):end,n);
nx = sqrt(nk/2);
x = reshape(x, nx, nx);
y = reshape(y, nx, nx);
p = y.^2+x.^2;
pm = p/max(p(:));

% spatial RF
imagesc(xax, yax, pm); colormap gray
hold on
contour(xax, yax, pm, .5, 'y');
axis xy
title('spatial STA')

subplot(1, 3, 2)
[xxx] = PlotSpatialk(v(:,n)*sflip, nx, spatres, fieldCenter, 'r'); hold on
PlotSpatialk(spatialK, nx, spatres, fieldCenter, 'k');
title('STA')

xlim(xxx(1,[1 end])+[-2 2])

subplot(1,3,3)
plot((0:(nkt-1))/frameRate, sflip*u(:,1), 'r', timeLags, temporalRF, 'k')
title('temporal STA')
legend({'STA', 'true'})


%% Plotting regression solution
% linear regression can overcome some stimulus correlations by whitening
% the input -- the covariance of the stimulus is inverted.
if ~exist('rho', 'var')
    rho = 2e1; % play with the ridge parameter
end
wmap = (xx + rho*eye(size(xx)))\xy;
W = flipud(reshape(wmap(1:end-addDC)/nY, nkt, nk));
[u,~,v] = svd(W);
n = 1;
sflip = sign(sum(u(1:floor(end/2),n)));

fig = figure(2); clf;
set(fig, 'Position', [1 515 1299 291]);

subplot(1, 3, 1)
x = v(1:end/2,n);
y = v(((end/2)+1):end,n);
nx = sqrt(nk/2);
x = reshape(x, nx, nx);
y = reshape(y, nx, nx);
p = y.^2+x.^2;
pm = p/max(p(:));

imagesc(xax, yax, pm); colormap gray
hold on
contour(xax, yax, pm, .5, 'y');
title('recovered spatial RF')

subplot(1, 3, 2)
[xxx] = PlotSpatialk(v(:,n)*sflip, nx, spatres, fieldCenter, 'r'); hold on
PlotSpatialk(spatialK, nx, spatres, fieldCenter, 'k');
title('recovered tuning')

xlim(xxx(1,[1 end])+[-2 2])

subplot(1,3,3)
plot((0:(nkt-1))/frameRate, sflip*u(:,1), 'r', timeLags, temporalRF, 'k')
suptitle(sprintf('Rho: %02.2f', rho))

%% PC regression
% principal components regression can overcome some correlations by
% projecting the stimulus onto the first m principal components and doing
% regression in the subspace



figure(38142); semilogy(diag(s), 'o-');
xlabel('\lambda')
ylabel('eigenvalue')


hold on
sd = diag(s);
semilogy(sd(1:nPrincipalComponents), 'or');
legend('throw-out', 'keep')



%%
% project back into velocity space
w = bas*wb;
% full PC regression solution
W = flipud(reshape(w(1:end-addDC)/nY, nkt, nk));
% space-time separable?
[u,~,v] = svd(W);
n = 1;
sflip = sign(sum(u(1:floor(end/2),n)));

fig = figure(2224); clf;
set(fig, 'Position', [1 515 1299 291]);

subplot(1, 3, 1)
x = v(1:end/2,n);
y = v(((end/2)+1):end,n);
nx = sqrt(nk/2);
x = reshape(x, nx, nx);
y = reshape(y, nx, nx);
p = y.^2+x.^2;
pm = p/max(p(:));

imagesc(xax, yax, pm); colormap gray
hold on
C = contour(xax, yax, pm, .5, 'y');

subplot(1, 3, 2)
[xxx] = PlotSpatialk(v(:,n)*sflip, nx, spatres, fieldCenter, 'r'); hold on
PlotSpatialk(spatialK, nx, spatres, fieldCenter, 'k');

xlim(xxx(1,[1 end])+[-2 2])

subplot(1,3,3)
plot((0:(nkt-1))/frameRate, sflip*u(:,1), 'r', timeLags, temporalRF, 'k')

suptitle(sprintf('# Principal Components used: %02.0f', nPrincipalComponents))