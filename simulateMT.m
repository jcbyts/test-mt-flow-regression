% parameters of our simulated cell
if ~exist('spatialSigma', 'var')
    spatialSigma = 2.5;
end

if ~exist('directionPreference', 'var')
    directionPreference        = 45;
end

if ~exist('nkt', 'var')
    nkt = 31;
end

% Load real data to test how this compares to what we show MT neurons
assert(exist('X', 'var')==1, 'you must have a stimulus created. Try testHyperflowField.m')

% Play movie of the stimulus flow field
figure(1); clf
for ii = 1:100
    PlotSpatialk(X(ii,:), ny);
    pause(.015)
end
    

%% make our receptive field (space time separable)
if ~exist('fieldCenter', 'var')
    fieldCenter = [0 0];
end

[xx,yy] = meshgrid(xax, yax);

if ~exist('RFcenter', 'var')
    RFcenter = fieldCenter + spatres/nx;
end

% make spatial RF a gaussian in velocity space
spatialRFx = cosd(directionPreference)*exp(- ((xx-RFcenter(1)).^2 + (yy-RFcenter(2)).^2).^.6/spatialSigma);
spatialRFy = sind(directionPreference)*exp(- ((xx-RFcenter(1)).^2 + (yy-RFcenter(2)).^2).^.6/spatialSigma);

% make time
timeLags = (0:30)/60;
tau = .05;
% alpha function... could be expanded to have a lag and a biphasic
% component
temporalRF = timeLags.*exp(-timeLags/tau)/tau;

% plot the simulated MT neuron
figure(1); clf
subplot(121)
spatialK = [spatialRFx(:); spatialRFy(:)];
PlotSpatialk(spatialK, size(xx));
xlabel('x space')
ylabel('y space')
title('simulated MT RF')
subplot(122)
plot(timeLags, temporalRF, 'k')
xlabel('time')
ylabel('amplitude')
title('temporal response MT')
suptitle('simulated "true" MT')

%% filter the stimulus with our simulated MT cell

% assume we have a space-time separable neuron
% project the velocity stimulus onto the spatial RF
rspace = X*spatialK;
% filter with the temporal RF
r  = filter(temporalRF, 1, rspace);

% pass through a nonlinearity to get spike rate
% we chose a compressive nonlinearity where firing rate saturates
nonlinearity = @(x, sc, in, st) sc./(1 + exp(-st*x+in));
nlinfun = @(x) nonlinearity(x, 4, 10, 4e-2);

xr = -10:.1:200; % range the nonlinearity takes

figure(2231); clf

subplot(131)
PlotSpatialk(spatialK, size(xx));
xlabel('x space')
ylabel('y space')
title('simulated MT RF')
subplot(132)
plot(xr, nlinfun(xr), 'k')
title('nonlinearity')
xlabel('k^Tx')
ylabel('rate')

if ~exist('mFr', 'var') % mean firing rate of simulated neuron
    mFr = 10;
end

% simulated responses
nl = nlinfun(r);
nl = nl/mean(nl) * mFr/60;

Y = poissrnd(nl);
subplot(133)
plot(nl./(1/60), 'k'); hold on
plot(find(Y), mFr*3, '.r')
xlim([1e3 1.5e3])   
xlabel('frame #')
ylabel('spike rate')
title('simulated neuron')
suptitle('simulated MT')
legend('rate', 'spikes')

clear r nl rspace
% %% Can we recover the neuron?
% 
% % build design matrix
% addDC = 1;
% nk = size(X,2);
% B = makeStimRows(X, nkt, 'same'); % builds the design matrix
% ny = numel(Y(:,1));
% B = [B ones(ny, addDC)];
% % pre compute som useful values
% xx = B'*B;
% yy = Y'*Y;
% xy = B'*Y;
% 
% %% plotting the Spike-Triggered Average
% % the spike triggered average is easy to compute and it is proportional to
% % the maximum-likelihood kernel if the stimulus is white noise. Ours isn't,
% % so you'll see that the spatial and temporal RFs seem blured out by the
% % correlations.
% 
% % this is the full spatiotemporal STA
% STA = flipud(reshape(xy(1:end-addDC)/ny, nkt, nk));
% 
% % check if it's 
% [u,~,v] = svd(STA);
% n = 1;
% sflip = sign(sum(u(1:floor(end/2),n)));
% 
% % plot the spike triggered average
% fig = figure(14924); clf;
% set(fig, 'Position', [1 515 1299 291]);
% subplot(1, 3, 1)
% x = v(1:end/2,n);
% y = v(((end/2)+1):end,n);
% nx = sqrt(nk/2);
% x = reshape(x, nx, nx);
% y = reshape(y, nx, nx);
% p = y.^2+x.^2;
% pm = p/max(p(:));
% 
% % spatial RF
% imagesc(xax, yax, pm); colormap gray
% hold on
% contour(xax, yax, pm, .5, 'y');
% axis xy
% title('spatial STA')
% 
% subplot(1, 3, 2)
% [xxx] = PlotSpatialk(v(:,n)*sflip, nx, spatres, fieldCenter, 'r'); hold on
% PlotSpatialk(spatialK, nx, spatres, fieldCenter, 'k');
% title('STA')
% 
% xlim(xxx(1,[1 end])+[-2 2])
% 
% subplot(1,3,3)
% plot(timeLags, sflip*u(:,1), 'r', timeLags, temporalRF, 'k')
% title('temporal STA')
% legend({'STA', 'true'})
% 
% 
% %% Plotting regression solution
% % linear regression can overcome some stimulus correlations by whitening
% % the input -- the covariance of the stimulus is inverted.
% if ~exist('rho', 'var')
%     rho = 2e1; % play with the ridge parameter
% end
% wmap = (xx + rho*eye(size(xx)))\xy;
% W = flipud(reshape(wmap(1:end-addDC)/ny, nkt, nk));
% [u,~,v] = svd(W);
% n = 1;
% sflip = sign(sum(u(1:floor(end/2),n)));
% 
% fig = figure(2); clf;
% set(fig, 'Position', [1 515 1299 291]);
% 
% subplot(1, 3, 1)
% x = v(1:end/2,n);
% y = v(((end/2)+1):end,n);
% nx = sqrt(nk/2);
% x = reshape(x, nx, nx);
% y = reshape(y, nx, nx);
% p = y.^2+x.^2;
% pm = p/max(p(:));
% 
% imagesc(xax, yax, pm); colormap gray
% hold on
% contour(xax, yax, pm, .5, 'y');
% title('recovered spatial RF')
% 
% subplot(1, 3, 2)
% [xxx] = PlotSpatialk(v(:,n)*sflip, nx, spatres, fieldCenter, 'r'); hold on
% PlotSpatialk(spatialK, nx, spatres, fieldCenter, 'k');
% title('recovered tuning')
% 
% xlim(xxx(1,[1 end])+[-2 2])
% 
% subplot(1,3,3)
% plot(timeLags, sflip*u(:,1), 'r', timeLags, temporalRF, 'k')
% 
% %% PC regression
% % principal components regression can overcome some correlations by
% % projecting the stimulus onto the first m principal components and doing
% % regression in the subspace
% 
% [b,s] = svd(xx);
% 
% figure(38142); semilogy(diag(s), 'o-');
% xlabel('\lambda')
% ylabel('eigenvalue')
% 
% if ~exist('nPrincipalComponents', 'var')
%     nPrincipalComponents = 1000;
% end
% hold on
% sd = diag(s);
% semilogy(sd(1:nPrincipalComponents), 'or');
% legend('throw-out', 'keep')
% 
% bas = b(:,1:nPrincipalComponents);
% 
% Xb = B*bas;
% Xb = bsxfun(@minus, Xb, mean(Xb));
% 
% xxb = Xb'*Xb;
% xyb = Xb'*Y;
% 
% wb  = xxb\xyb;
% %%
% % project back into velocity space
% w = bas*wb;
% % full PC regression solution
% W = flipud(reshape(w(1:end-addDC)/ny, nkt, nk));
% % space-time separable?
% [u,~,v] = svd(W);
% n = 1;
% sflip = sign(sum(u(1:floor(end/2),n)));
% 
% fig = figure(2); clf;
% set(fig, 'Position', [1 515 1299 291]);
% 
% subplot(1, 3, 1)
% x = v(1:end/2,n);
% y = v(((end/2)+1):end,n);
% nx = sqrt(nk/2);
% x = reshape(x, nx, nx);
% y = reshape(y, nx, nx);
% p = y.^2+x.^2;
% pm = p/max(p(:));
% 
% imagesc(xax, yax, pm); colormap gray
% hold on
% C = contour(xax, yax, pm, .5, 'y');
% 
% subplot(1, 3, 2)
% [xxx] = PlotSpatialk(v(:,n)*sflip, nx, spatres, fieldCenter, 'r'); hold on
% PlotSpatialk(spatialK, nx, spatres, fieldCenter, 'k');
% 
% xlim(xxx(1,[1 end])+[-2 2])
% 
% subplot(1,3,3)
% plot(timeLags, sflip*u(:,1), 'r', timeLags, temporalRF, 'k')

% %% use cosine basis for time
% 
% % The STA here is working out well... we have to figure out what to do with
% % the regression bit though. It's not working :(
% nb = 13;
% nk = size(X,2);
% kUnit = 5;
% addDC = 1;
% [iht, ihbas, ihbasis] = makeRaisedCosBasis(nb, 1, [0 nkt], 10);
% 
% 
% figure(2); clf
% plot(iht, ihbasis);
% 
% tbx = @(x) temporalBases_dense(x, ihbasis, ones(nk, nb), addDC);
% 
% 
% zaso = encapsulateRaw(X, Y, tbx);
% 
% fnlin = @expfun;
% 
% [rsum, ragg] = zasoFarray(zaso, {@(x,y) x' * x, @(x,y) x' * y, @(x,y) y'*y}, {});
% 
% 
% %%
% map = rsum{2};
% map = (rsum{1} + 8e5*eye(size(rsum{1})))\rsum{2};
% A = reshape(map(1:end-addDC), nb, []);
% 
% At = ihbasis*A;
% 
% figure(13454); clf
% imagesc(At)
% 
% [u,~,v] = svd(At);
% n = 1;
% sflip = sign(sum(u(1:floor(end/2),n)));
% 
% fig = figure(2); clf;
% set(fig, 'Position', [1 515 1299 291]);
% 
% subplot(1, 3, 1)
% x = v(1:end/2,n);
% y = v(((end/2)+1):end,n);
% nx = sqrt(nk/2);
% x = reshape(x, nx, nx);
% y = reshape(y, nx, nx);
% p = y.^2+x.^2;
% pm = p/max(p(:));
% 
% imagesc(xax, yax, pm); colormap gray
% hold on
% C = contour(xax, yax, pm, .5, 'y');
% 
% subplot(1, 3, 2)
% [xxx] = hyperflow.PlotSpatialk(v(:,n)*sflip, nx, spatres, fieldCenter, 'r'); hold on
% hyperflow.PlotSpatialk(spatialK, nx, spatres, fieldCenter, 'k');
% 
% xlim(xxx(1,[1 end])+[-2 2])
% 
% subplot(1,3,3)
% plot(iht/60, sflip*u(:,1), 'r', timeLags, temporalRF, 'k')

