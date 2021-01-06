%{
------------HEADER-----------------
Objective          ::  Estimate the noise covariance matrix
INPUT VARS
trs              :: The telemetry top-level class
nMin, nMax  :: (optional) 'Min and Max Zernike modes order the wavefront is decomposed onto (defaut: 4, 100)
fitL0             :: (optional) fit the L0 as well if set to true (default: true)
best              :: (optional) Repeat the fitting process for different values of nMin from nMin to nMax/4 and pick up the most precise one
wvl_ref         :: (optional) Wavelength the r0 is given at.
mskPhase     :: (optional) Phase mask to reject non-valid actuators.
OUTPUT VARS
r0               :: The estimated line-of-sight r0 in meters at wvl
L0               :: The estimated L0 in meters
fwhm          :: The estimated PSF FWHM in arcsec
dr0             :: The 3-sigma precision on the r0 estimation
dL0             :: The 3-sigma precision on the L0 estimation
dfwhm        :: The 3-sigma precision on the fwhm estimation
Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}


function [resSeeing,resZer] = estimateSeeingFromTelemetry(trs,varargin)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));
inputs.addParameter('wvl',0.5e-6,@isnumeric);
inputs.addParameter('fitL0',true,@islogical);
inputs.addParameter('flagBest',false,@islogical);
inputs.addParameter('flagMedian',false,@islogical);
inputs.addParameter('D1',9,@isnumeric);
inputs.addParameter('D2',2.65,@isnumeric);
inputs.addParameter('badModesList',[],@isnumeric);
inputs.addParameter('aoMode','NGS',@ischar);
inputs.addParameter('jMin',4,@isnumeric);
inputs.addParameter('jMax',120,@isnumeric);
inputs.parse(trs,varargin{:});


%1\ Parsing inputs 
fitL0 = inputs.Results.fitL0;
flagBest = inputs.Results.flagBest;
flagMedian = inputs.Results.flagMedian;
D1 = inputs.Results.D1; % 9 Equivalent outer diameter of a circular pupil within which actuators commands are not sensitive to edges effects
D2 = inputs.Results.D2; % 2.65 Same for the inner diameter
badModesList = inputs.Results.badModesList;
aoMode = inputs.Results.aoMode;
jMin = inputs.Results.jMin;
jMax = inputs.Results.jMax;
wvl = inputs.Results.wvl;
if isfield(trs.mat,'dmIF_hr') && ~isempty(trs.mat.dmIF_hr)
    dmModes = trs.mat.dmIF_hr;
else
    dmModes = trs.mat.dmIF;
end

%2\ Getting the reconstructed wavefront and noise covariance at 500 nm
[uout,mskModes] = selectActuators(trs.dm.com,D1,D2);
uout = (2*pi/wvl)*uout;
validActu = mskModes(:);
nValid = nnz(mskModes);

%3\ Getting the pupil mask
nRes = sqrt(size(dmModes,1));
res = getGridCoordinates(nRes,nRes,10.54*trs.dm.pitch);
mskOuter = abs(res.x2D) <= D1/2 & abs(res.y2D) <= D1/2;
mskInner = res.r2D <= D1/2 & res.r2D >= D2/2;
dmModes = bsxfun(@times,dmModes,mskInner(:));


%3\ Get the noise covariance matrix and variance
if ~isfield(trs.res,'noise') || (isfield(trs.res,'noise')  && isempty(trs.res.noise(1).Cn_ho))
    fprintf('Estimate the noise covariance matrix\n')
    getNoiseSTD(trs);
end
Cn   = 0*trs.res.noise(1).Cn_ho;
Cn2 = trs.res.noise(1).Cn_ho;
Cn(validActu,validActu) = Cn2(validActu,validActu)*(2*pi/wvl)^2;
varN = trace(Cn(validActu,validActu))/nValid;

%4\ Preliminary r0 estimation sig^2 = 0.111*(D(1-o)/r0)^5/3
varPh = sum(std(uout(validActu,:),[],2).^2)/nValid - varN;
jMaxData = round((sqrt(nValid)+1)*(sqrt(nValid)+2)/2-1);
kAO = 0.134;%piston and tip-tilt excluded
kFit = 0.2944*jMaxData^(-sqrt(3)/2);
%var = (kAO - kFit)*(D/r0)^5/3
r0_ = min(D1/(varPh/(kAO-kFit))^(3/5),0.2);


%4\ Estimating the r0/outer scale
opt = {'fitL0',fitL0,'flagBest',flagBest,'flagMedian',flagMedian,'initR0',r0_,'D1',D1,'D2',D2,...
    'badModesList',badModesList,'aoMode',aoMode,'jMin',jMin,'jMax',jMax,'mskModes',mskModes,'mskPup',mskOuter};
[r0,L0,dr0,dL0,jind,z_mod,z_meas,z_noise] =  getr0L0FromDMcommands(dmModes,uout,trs.tel,Cn,opt{:});

resZer.fitL0 = fitL0;
resZer.flagBest = flagBest;
resZer.flagMedian = flagMedian;
resZer.r0_init = r0_;


% !!!! IF the precision is not good,adjust the r0 only !!!
if dr0*3 > 1
    opt = {'fitL0',false,'flagBest',flagBest,'flagMedian',flagMedian,'initR0',r0_,'D1',D1,'D2',D2,...
        'badModesList',badModesList,'aoMode',aoMode,'jMin',jMin,'jMax',jMax,'mskModes',mskModes,'mskPup',mskOuter};
    [r0,~,dr0,~,jind,z_mod,z_meas,z_noise] =  getr0L0FromDMcommands(dmModes,uout,trs.tel,Cn,opt{:});
end

fprintf('r_0 estimated from the telemetry: %.3g cm\n',1e2*r0);
fprintf('L_0 estimated from the telemetry: %.3g m\n',L0);

%5\ Estimating the atmosphere-limited PSF FWHM
if r0~=r0
    r0 = r0_;
end
k1 = 0.976*wvl*constants.radian2arcsec;
seeing = k1/r0;
dseeing= k1*dr0/r0^2;

if L0 ~=0
    k2     = 2.813;
    k3     = 0.356;
    sqfact = sqrt(1-k2*(r0/L0)^k3);
    seeing = seeing*sqfact;
    
    dfdr01 = dr0/r0*seeing;
    dfdr02 = 0.5*k1*dr0/r0/sqfact*k2*k3*L0^(-k3)*r0^(k3-1);
    dfdr0 = hypot(dfdr01,dfdr02);        
    dfdL0 = 0.5*k1*dL0/r0/sqfact*k2*k3*r0^(k3)*L0^(-k3-1);    
    dseeing=hypot(dfdr0,dfdL0);
end

% Concatenating in structure results
resSeeing.r0 = r0;
resSeeing.L0 = L0;
resSeeing.seeing_VonKarman = seeing;
resSeeing.dr0 = dr0;
resSeeing.dL0 = dL0;
resSeeing.dseeing_VonKarman = dseeing;
resSeeing.seeing_Kolmo = 0.976*constants.radian2arcsec*0.5e-6/r0;
resSeeing.dseeing_Kolmo = dr0*0.976*constants.radian2arcsec*0.5e-6/r0^2;
resZer.jindex = jind;
resZer.std_model = z_mod;
resZer.std_meas = z_meas;
resZer.std_noise = z_noise;


