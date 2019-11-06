function psfr = forwardPSFR(psfr,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfReconstruction'));
inputs.addParameter('flagToeplitz',false,@islogical);
inputs.addParameter('flagAnisoMethod','oomao',@ischar);
inputs.addParameter('flagNoiseMethod','autocorrelation',@ischar);
inputs.addParameter('flagAoPattern','circle',@ischar);
inputs.parse(psfr,varargin{:});

psfr.flags.toeplitz         = inputs.Results.flagToeplitz;
psfr.flags.anisoMethod = inputs.Results.flagAnisoMethod;
psfr.flags.noiseMethod = inputs.Results.flagNoiseMethod;
psfr.flags.aoPattern = inputs.Results.flagAoPattern;

%1\ Noise covariance matrices estimation
if (isprop(psfr.trs,'res') || isfield(psfr.trs,'res')) && isfield(psfr.trs.res,'noise') && ~isempty(psfr.trs.res.noise)  && (~isempty(psfr.trs.res.noise.Cn_ho) && ~isempty(psfr.trs.res.noise.Cn_tt))
    fprintf('The noise covariance matrix has been already successfully estimated.\n');
else
    fprintf('Estimating the noise covariance matrix\n');
    psfr.trs.res.noise = estimateNoiseCovarianceFromTelemetry(psfr.trs,'method',psfr.flags.noiseMethod);
end

%\2 Seeing estimation
if (isprop(psfr.trs,'res') || isfield(psfr.trs,'res')) && isfield(psfr.trs.res,'seeing') && ~isempty(psfr.trs.res.seeing) && ~isempty(psfr.trs.res.seeing.r0) 
    fprintf('The seeing has been already successfully estimated.\n');
else
    fprintf('Estimating the seeing value\n');
    [psfr.trs.res.seeing,psfr.trs.res.zernike] = estimateSeeingFromTelemetry(psfr.trs);
end

%3\ Diffraction-limit OTF - Nyquist sampling
psfr.otf.otfStat= computeStaticOpticalTransferFunction(psfr);

%4\ Normalized Fitting SF
[psfr.sf.Dfit,psfr.cov.psdFit] = computeFittingPhaseStructureFunction(psfr,'aoPattern',psfr.flags.aoPattern);

%5\ Normalized Alasing SF
[psfr.sf.Dal,psfr.cov.Cal] = computeAliasingPhaseStructureFunction(psfr,'aoPattern',psfr.flags.aoPattern);

%6\ AO residual SF
[psfr.sf.Dho_z,psfr.cov.Cho_z] = computeResidualPhaseStructureFunction(psfr);
psfr.sf.Dho = psfr.sf.Dho_z;
psfr.cov.Cho = psfr.cov.Cho_z;

%7\ Tip-tilt SF
[psfr.sf.Dtt,psfr.cov.Ctt] = computeTipTiltPhaseStructureFunction(psfr);

%8\ Anisoplanatism
[psfr.sf.Dani_l,psfr.sf.Dani,psfr.otf.Kani] = computeAnisoplanatismPhaseStructureFunction(psfr);

%9\  transfer function
psfr.otf.otfCCD = computeCcdOpticalTransferFunction(psfr.trs.cam.pixelScale*1e-3,2*psfr.trs.tel.Dcircle/psfr.trs.cam.wavelength/psfr.otf.nOtf,psfr.otf.nOtf);

%9\ Reconstruct the PSF
psfr.psf.rec = zeros(psfr.psf.fov,psfr.psf.fov,psfr.trs.src.nSrc);
for iSrc = 1:psfr.trs.src.nSrc
    % Get the OTFwith a PSF Nyquist-sampling
    psfr.otf.otfShannon = psfr.otf.otfStat.*exp(-0.5*(psfr.sf.Dfit+ psfr.sf.Dal+ psfr.sf.Dho + psfr.sf.Dtt + psfr.sf.Dani(:,:,iSrc))).*psfr.otf.otfCCD;  
    % Get the PSF at the imager pixel scale
    if abs(psfr.trs.tel.pupilAngle) > 2
        psfr.otf.otfShannon = rotateImage(psfr.otf.otfShannon,psfr.trs.tel.pupilAngle);
    end
    psf_ij = tools.otfShannon2psf(psfr.otf.otfShannon,psfr.trs.cam.samp(iSrc),psfr.psf.fov);
    % Normalization across the support the PSF will be cropped    
    S = sum(sum(tools.crop(psf_ij,psfr.trs.cam.resolution)));
    psfr.psf.rec(:,:,iSrc) = psf_ij/S;    
end

