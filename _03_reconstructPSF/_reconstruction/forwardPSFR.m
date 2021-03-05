function psfr = forwardPSFR(psfr,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfReconstruction'));
inputs.addParameter('flagToeplitz',false,@islogical);
inputs.addParameter('flagAoPattern','circle',@ischar);
inputs.addParameter('flagAnisoMethod','flicker',@ischar);
inputs.addParameter('flagNoiseMethod','autocorrelation',@ischar);
inputs.addParameter('flagDphiMethod','Vii',@ischar);
inputs.addParameter('flagResidualMethod','slopes-based',@ischar);
inputs.addParameter('statModes',{[]},@iscell);
inputs.parse(psfr,varargin{:});

psfr.flags.toeplitz      = inputs.Results.flagToeplitz;
psfr.flags.anisoMethod   = inputs.Results.flagAnisoMethod;
psfr.flags.noiseMethod   = inputs.Results.flagNoiseMethod;
psfr.flags.aoPattern     = inputs.Results.flagAoPattern;
psfr.flags.dphiMethod    = inputs.Results.flagDphiMethod;
psfr.flags.residualMethod= inputs.Results.flagResidualMethod;
psfr.statModes           = inputs.Results.statModes;

%1\ Noise covariance matrices estimation
if (isprop(psfr.trs,'res') || isfield(psfr.trs,'res')) && isfield(psfr.trs.res,'noise') && ~isempty(psfr.trs.res.noise) ...
        && (~isempty(psfr.trs.res.noise.Cn_ho) && ~isempty(psfr.trs.res.noise.Cn_tt)) && ...
        strcmp(psfr.flags.noiseMethod,psfr.trs.res.noise.method)
    fprintf('The noise covariance matrix has been already successfully estimated.\n');
else
    fprintf('Estimating the noise covariance matrix\n');
    psfr.trs.res.noise = estimateNoiseCovarianceFromTelemetry(psfr.trs,'flagNoisemethod',psfr.flags.noiseMethod);
end

%\2 Seeing estimation
if (isprop(psfr.trs,'res') || isfield(psfr.trs,'res')) && isfield(psfr.trs.res,'seeing') && ~isempty(psfr.trs.res.seeing) && ~isempty(psfr.trs.res.seeing.r0) 
    fprintf('The seeing has been already successfully estimated.\n');
else
    fprintf('Estimating the seeing value\n');
    [psfr.trs.res.seeing,psfr.trs.res.zernike] = estimateSeeingFromTelemetry(psfr.trs);
end

%3\ Diffraction-limit OTF - Nyquist sampling
psfr = computeStaticOpticalTransferFunction(psfr,'statModes',psfr.statModes);

%4\ Normalized Fitting SF
psfr = computeFittingPhaseStructureFunction(psfr,'aoPattern',psfr.flags.aoPattern);

%5\ Normalized Alasing SF
psfr = computeAliasingPhaseStructureFunction(psfr,'aoPattern',psfr.flags.aoPattern);

%6\ AO residual SF
psfr = computeResidualPhaseStructureFunction(psfr,'method',psfr.flags.dphiMethod,'flagResidualMethod',psfr.flags.residualMethod);

%7\ Tip-tilt SF
psfr = computeTipTiltPhaseStructureFunction(psfr);

%8\ Anisoplanatism
psfr = computeAnisoplanatismPhaseStructureFunction(psfr);

%9\ CCD  transfer function 
if ~isempty(psfr.trs.obj_name)
    psfr.otf.otfCCD = computeCcdOpticalTransferFunction(psfr.trs.cam.pixelScale*1e-3,2*psfr.trs.tel.Dcircle/psfr.trs.cam.wavelength/psfr.otf.nOtf,psfr.otf.nOtf);
else
    %Simulation case
    psfr.otf.otfCCD = 1;
end

%9\ Reconstruct the PSF
psfr.rec_ = zeros(psfr.otf.fov_fit,psfr.otf.fov_fit,numel(psfr.trs.src));


for iSrc = 1:numel(psfr.trs.src)
    % Get the OTFwith a PSF Nyquist-sampling
    if contains(psfr.flags.dphiMethod,'zonal')
        psfr.otf.otfShannon = psfr.otf.otfStat.*psfr.otf.otfAO.*exp(-0.5*(psfr.sf.Dfit+ psfr.sf.Dal+ psfr.sf.Dtt + psfr.sf.Dani(:,:,iSrc))).*psfr.otf.otfCCD;
    else
        psfr.otf.otfShannon = psfr.otf.otfStat.*exp(-0.5*(psfr.sf.Dfit+ psfr.sf.Dal+ psfr.sf.Dao + psfr.sf.Dtt + psfr.sf.Dani(:,:,iSrc))).*psfr.otf.otfCCD;
    end
    % Get the PSF at the imager pixel scale   
    psf_ij = puakoTools.otfShannon2psf(psfr.otf.otfShannon,psfr.trs.cam.samp,psfr.otf.fov_fit);     
    % Normalization across the support the PSF will be cropped    
    S = sum(sum(puakoTools.crop(psf_ij,psfr.trs.cam.resolution)));
    psfr.rec_(:,:,iSrc) = psf_ij/S;        
end

