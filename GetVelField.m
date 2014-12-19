function stimfull = GetVelField(params, opticflow, aperturecenter)
% Reconstruct Velocity Field from optic flow seeds
% Input:
%       params: parameters for the velocity field
%               - designsizex designsizey, size of velocity field (in unit
%               of degrees, e.g., 30)
%               - spatres, spatial resolution (in units of degrees, e.g. 2)
%               dimensions
%               - maskdiameter, diameter of the aperture
%       opticflow: NTx6 six optic flow components
%       aperturecenter: NTx2 aperture center trajectory
% Yuwei Cui Dec 2, 2013
newdatalen = size(opticflow,1);

params.samplingratex = 1/params.spatres;
params.samplingratey = 1/params.spatres;

xl = round(params.designsizex/params.spatres);
yl = round(params.designsizey/params.spatres);

[xi0, yi0] = meshgrid(((-xl/2+0.5):(xl/2-0.5))/params.samplingratex,...
    ((-yl/2+0.5):(yl/2-0.5))/params.samplingratey);

% aperture center
centerx = aperturecenter(:,1);
centery = aperturecenter(:,2);

stimfull = zeros(newdatalen,xl*yl*2);
for ii = 1:newdatalen
    xi = xi0-centerx(ii);
    yi = yi0-centery(ii);
    % make mask
    mask = (xi.^2+yi.^2<=(params.maskdiameter/2)^2);
    
    rawbasismat = zeros(round(params.designsizey*params.samplingratey),...
        round(params.designsizex*params.samplingratex),2,6);
    % translation
    rawbasismat(:,:,1,1) = 1.*mask;
    rawbasismat(:,:,2,1) = 0;
    rawbasismat(:,:,1,2) = 0;
    rawbasismat(:,:,2,2) = 1.*mask;
    % expansion
    rawbasismat(:,:,1,3) = xi.*mask;
    rawbasismat(:,:,2,3) = yi.*mask;
    % rotation
    rawbasismat(:,:,1,4) = -yi.*mask;
    rawbasismat(:,:,2,4) = xi.*mask;
    % shear 1
    rawbasismat(:,:,1,5) = xi.*mask;
    rawbasismat(:,:,2,5) = -yi.*mask;
    % shear 2
    rawbasismat(:,:,1,6) = yi.*mask;
    rawbasismat(:,:,2,6) = xi.*mask;
    %                                     keyboard
    % reconstruct stimuli
    stimfull(ii,:) = reshape(rawbasismat,[],6)*opticflow(ii,1:6)';
    
end
stimfull(isnan(stimfull))=0;


% % flip stimulus
% stimflip = zeros(size(stimfull));
% Nx = sqrt(size(stimflip,2)/2); Ny = Nx;
% for i=1:size(stimflip,1)
%     VelField = stimfull(i,:);
%     Vx = VelField(1:end/2); Vy = VelField(end/2+1:end);
%     Vx = flipud(reshape(Vx, Nx, Ny)); Vy = flipud(reshape(Vy, Nx, Ny));
%     stimflip(i,:) = [Vx(:); Vy(:)];
% end
% stimfull = stimflip;
% end