function trs = fromAoSystemClassToTelemetry(trs,aoSys)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));
inputs.addRequired('aoSys',@(x) isa(x,'aoSystem'));

rad2arcsec = 3600 * 180/pi;
%% 1\ Restore simulations settings
%1.1 AO config
trs.aoMode = 'NGS';
if isprop(aoSys,'lGs')
    trs.aoMode = 'LGS';    
end

%1.2. Atmospheric config
trs.atm.wavelength      = aoSys.atm.wavelength;
trs.atm.r0              = aoSys.atm.r0;
trs.atm.L0              = aoSys.atm.L0;
trs.atm.weights         = [aoSys.atm.layer.fractionnalR0];
trs.atm.heights         = [aoSys.atm.layer.altitude];
trs.atm.Cn2             = trs.atm.r0 ^(-5/3)*trs.atm.weights;
trs.atm.windSpeed       = [aoSys.atm.layer.windSpeed];
trs.atm.windDirection   = [aoSys.atm.layer.windDirection];
trs.atm.nLayer          = length(trs.atm.weights);
trs.atm.zenith_angle    = aoSys.parms.atm.zenithAngle;

%1.3. Pupil
trs.tel.Ddm             = aoSys.tel.D;
trs.tel.pupil           = aoSys.tel.pupil;
trs.tel.resolution      = aoSys.tel.resolution;
trs.tel.pixelscale      = aoSys.tel.D/trs.tel.resolution;
trs.tel.airmass         = 1/cos(trs.atm.zenith_angle*pi/180);
trs.tel.pupilAngle      = 0;

%1.4. Static aberrations
trs.cam.name            = 'OOMAO';
trs.tel.static_map      = aoSys.ncpa.map;
trs.cam.x_ncpa          = 0;
trs.cam.y_ncpa          = 0;
trs.tel.static_fitting  = [];

%% 2\ Restore AO telemetry
%2.1. WFS slopes
trs.wfs.slopes  = aoSys.matrices.SlopeTTRem*squeeze(aoSys.loopData.slopes);
trs.wfs.slopes  = bsxfun(@minus,trs.wfs.slopes ,mean(trs.wfs.slopes ,2));
trs.wfs.nSl     = size(trs.wfs.slopes,1);
trs.wfs.nExp    = size(trs.wfs.slopes,2);
trs.wfs.theta   = 0;

trs.wfs.pixelScale = 0;
if isprop(aoSys,'lGs')
    trs.wfs.wavelength = aoSys.lGs.wavelength;
    trs.tipTilt.wavelength = aoSys.nGs.wavelength;
   
    % LGS source
    trs.lgs.height      = aoSys.lGs.height;
    xLGS                = [aoSys.lGs.directionVector];
    trs.lgs.x           = xLGS(1,:);
    trs.lgs.y           = xLGS(2,:);
    % LGS WFS pixel scale calculation
    d                   = aoSys.tel.D/size(aoSys.wfs.validLenslet,1);
    nPxDetector         = size(aoSys.wfs.camera.frame,1)/aoSys.wfs.lenslets.nLenslet;
    binFactor           = 2*aoSys.wfs.lenslets.fieldStopSize*aoSys.wfs.lenslets.nyquistSampling/nPxDetector;
    lo2DInMas           = aoSys.lGs(1).wavelength/(2*d)*constants.radian2mas;
    trs.wfs.pixelScale  = lo2DInMas/aoSys.wfs.lenslets.nyquistSampling*binFactor;
    
    % LO WFS pixel scale calculation : 
    d                   = aoSys.tel.D/size(aoSys.lowfs.validLenslet,1);
    nPxDetector         = size(aoSys.lowfs.camera.frame,1)/aoSys.lowfs.lenslets.nLenslet;
    binFactor           = 2*aoSys.lowfs.lenslets.fieldStopSize*aoSys.lowfs.lenslets.nyquistSampling/nPxDetector;
    lo2DInMas           = aoSys.nGs(1).wavelength/(2*d)*constants.radian2mas;
    trs.wfs.pixelScale  = lo2DInMas/aoSys.lowfs.lenslets.nyquistSampling*binFactor;
    
else
    trs.wfs.wavelength  = aoSys.nGs.wavelength;
    % NGS WFS pixel scale calculation
    d                   = aoSys.tel.D/size(aoSys.wfs.validLenslet,1);
    nPxDetector         = size(aoSys.wfs.camera.frame,1)/aoSys.wfs.lenslets.nLenslet;
    binFactor           = 2*aoSys.wfs.lenslets.fieldStopSize*aoSys.wfs.lenslets.nyquistSampling/nPxDetector;
    lo2DInMas           = aoSys.nGs(1).wavelength/(2*d)*constants.radian2mas;
    trs.wfs.pixelScale  = lo2DInMas/aoSys.wfs.lenslets.nyquistSampling*binFactor;
end

% NGS source
xNGS      = [aoSys.nGs.directionVector];
trs.ngs.x = xNGS(1,:);
trs.ngs.y = xNGS(2,:);
    
%2.2 WFS pixels intensity
nPix = aoSys.wfs.camera.resolution(1);
trs.wfs.intensity = zeros(nPix,nPix);
if isprop(aoSys.loopData,'intensity')
    trs.wfs.intensity = median(aoSys.loopData.intensity,3);
end

%2.3. DM commands
trs.dm.volt2meter       = 2;
trs.dm.com              = aoSys.matrices.DMTTRem*squeeze(aoSys.loopData.dmcom)*trs.dm.volt2meter;
trs.dm.com              = bsxfun(@minus,trs.dm.com,mean(trs.dm.com,2));
trs.dm.nCom             = size(trs.dm.com,1);
trs.dm.validActuators   = aoSys.dm.validActuator;
trs.dm.nActuators       = size(trs.dm.validActuators,1);
trs.dm.pitch            = aoSys.parms.dm.pitch;

%2.4 Influence DM functions
trs.dm.pupilMask    = logical(puakoTools.interpolate(trs.tel.pupil,trs.dm.nActuators));
bif                 = xineticsInfluenceFunctionModel(trs.dm.pitch,'mechCoupling',aoSys.dm.modes.mechCoupling);
bif.setInfluenceFunction(trs.dm.nActuators,2*trs.dm.nActuators-1);
trs.mat.dmIF        = bif.modes;
trs.mat.dmIF_inv    = pinv(full(trs.mat.dmIF));
trs.mat.Hdm         = trs.mat.dmIF*trs.mat.dmIF_inv;

%2.5. Tip-tilt
trs.tipTilt.tilt2meter  = 1;%!!!!!!
trs.tipTilt.slopes  = trs.tipTilt.tilt2meter*aoSys.loopData.tiptilt;
trs.tipTilt.slopes  = bsxfun(@minus,trs.tipTilt.slopes,mean(trs.tipTilt.slopes,3));

trs.tipTilt.com     = trs.tipTilt.tilt2meter*aoSys.loopData.tiltCom;
trs.tipTilt.com     = bsxfun(@minus,trs.tipTilt.com,mean(trs.tipTilt.com,2));

if isprop(aoSys,'lGs')
    trs.tipTilt.nExp    = size(trs.tipTilt.slopes,3);
else
    trs.tipTilt.nExp    = size(trs.tipTilt.slopes,2);
end
        
%% 3\ Get system matrices
%3.1\ Get DM commands reconstructors from slopes
MC          = aoSys.matrices.DMTTRem*aoSys.wfs2dm.M*aoSys.matrices.SlopeTTRem; %command matrix
trs.mat.R   = trs.dm.volt2meter*MC;
trs.mat.Rtt = trs.tipTilt.tilt2meter;

%3.2\ Get the reconstructed wavefront in OPD and in the actuators space
trs.rec.res = trs.dm.volt2meter*trs.mat.R*trs.wfs.slopes;
trs.rec.res = bsxfun(@minus,trs.rec.res,mean(trs.rec.res,2));

%3.3 fill vector to get 21x21 actuators
u                           = zeros(trs.dm.nActuators^2,trs.wfs.nExp);
u(trs.dm.validActuators,:)  = trs.rec.res;
trs.rec.res                 = u;
u                           = zeros(trs.dm.nActuators^2,trs.wfs.nExp);
u(trs.dm.validActuators,:)  = aoSys.matrices.DMTTRem*trs.dm.com;
trs.dm.com                  = u;

%% 4\ Get the loop status and model transfer function
%4.1. Latency
if ~isfield(aoSys.parms.loopStatus,'tt')
    aoSys.loopStatus.tt = aoSys.loopStatus.ho;
end
trs.holoop.lat      = aoSys.loopStatus.ho.latency;
trs.ttloop.lat      = aoSys.loopStatus.tt.latency;

%4.2. Frequency
trs.holoop.freq     = 1/aoSys.loopStatus.ho.frameRate/aoSys.tel.samplingTime;
trs.ttloop.freq     = 1/aoSys.loopStatus.tt.frameRate/aoSys.tel.samplingTime;

%4.3. RTC controller HO loop: pure proportional-integrator
trs.holoop.gain     = aoSys.loopStatus.ho.gain;
trs.holoop.tf.num   = [trs.holoop.gain 0 0 0];
trs.holoop.tf.den   = [-1 0 0];

trs.ttloop.gain     = aoSys.loopStatus.tt.gain;
trs.ttloop.tf.num   = [trs.ttloop.gain 0 0 0];
trs.ttloop.tf.den   = [-1 0 0];

%% 5\ Restore the PSF
tmp                 = aoSys.cam.imgLens.pixelScale(aoSys.sci(1),aoSys.tel); %recalculate the true pixel scale
trs.cam.pixelScale  = tmp.MAS;
trs.cam.resolution  = size(aoSys.psf.image,1);
trs.cam.fov         = trs.cam.resolution*trs.cam.pixelScale/1e3;

trs.cam.image       = aoSys.psf.image ;
trs.cam.wavelength  = aoSys.sci.wavelength;
trs.cam.samp        = constants.radian2mas*trs.cam.wavelength/trs.tel.Ddm/2/trs.cam.pixelScale;
trs.cam.itime       = aoSys.cam.exposureTime*aoSys.tel.samplingTime;
trs.cam.coadds      = 1;

% this assumes there is a single PSF at the center of the image
[trs.src.x,trs.src.y] = pol2cart(aoSys.sci.azimuth,aoSys.sci.zenith*rad2arcsec);
trs.src.F            = 10^(-0.4 *aoSys.sci.magnitude) * aoSys.sci.photometry.zeroPoint;

trs.sky = psfStats( trs.cam.image,trs.tel.pupil,trs.cam.wavelength,...
    trs.cam.samp,trs.cam.pixelScale,trs.src,'flagMoffat',false,'flagGaussian',false);
end