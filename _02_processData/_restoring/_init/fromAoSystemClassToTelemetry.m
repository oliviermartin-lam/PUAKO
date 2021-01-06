function trs = fromAoSystemClassToTelemetry(trs,aoSys)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));
inputs.addRequired('aoSys',@(x) isa(x,'aoSystem'));

%% 1\ Restore simulations settings
%1.1 AO config
trs.aoMode = 'NGS';
if ~isempty(aoSys.lGs)
    trs.aoMode = 'LGS';    
end
%1.2. Atmospheric config


trs.atm.wavelength = aoSys.atm.wavelength;
trs.atm.r0 = aoSys.atm.r0;
trs.atm.L0 = 25;
trs.atm.weights = [aoSys.atm.layer.fractionnalR0];
trs.atm.heights = [aoSys.atm.layer.altitude];
trs.atm.Cn2 = trs.atm.r0 ^(-5/3)*trs.atm.weights;
trs.atm.windSpeed =  [aoSys.atm.layer.windSpeed];
trs.atm.windDirection = [aoSys.atm.layer.windDirection];
trs.atm.nLayer = length(trs.atm.weights);

%1.3. Pupil
trs.tel.Ddm = aoSys.tel.D;
trs.tel.pupil = aoSys.tel.pupil;
trs.tel.resolution = aoSys.tel.resolution;
trs.tel.pixelscale = aoSys.tel.D/trs.tel.resolution;
trs.tel.zenith_angle = aoSys.parm.atm.zenithAngle;
trs.tel.airmass = 1/cos(trs.tel.zenith_angle*pi/180);
trs.tel.pupilAngle  = 0;

%1.4. Static aberrations
trs.cam.name = 'OOMAO';
trs.tel.static_map = aoSys.ncpa.map;
trs.cam.x_ncpa = 0;
trs.cam.y_ncpa = 0;
trs.tel.static_fitting = [];

%% 2\ Restore AO telemetry
%2.1. WFS slopes
trs.wfs.slopes = aoSys.matrices.SlopeTTRem*squeeze(aoSys.loopData.slopes);
trs.wfs.slopes     = bsxfun(@minus,trs.wfs.slopes ,mean(trs.wfs.slopes ,2));
trs.wfs.nSl     = size(trs.wfs.slopes,1);
trs.wfs.nExp   = size(trs.wfs.slopes,2);
trs.wfs.theta = 0;

%2.2 WFS pixels intensity
nPix = aoSys.wfs.camera.resolution(1);
trs.wfs.intensity = zeros(nPix,nPix);
trs.wfs.intensity = median(aoSys.loopData.intensity,3);

%2.3. DM commands
trs.dm.volt2meter = 2;
trs.dm.com    = aoSys.matrices.DMTTRem*squeeze(aoSys.loopData.dmcom)*trs.dm.volt2meter;
trs.dm.com= bsxfun(@minus,trs.dm.com,mean(trs.dm.com,2));
trs.dm.nCom = size(trs.dm.com,1);
trs.dm.validActuators = aoSys.dm.validActuator;
trs.dm.nActuators = size(trs.dm.validActuators,1);
trs.dm.pitch = aoSys.parm.dm.pitch;

%2.4 Influence DM functions
trs.dm.pupilMask = logical(puakoTools.interpolate(trs.tel.pupil,trs.dm.nActuators));
bif      = xineticsInfluenceFunctionModel(trs.dm.pitch,'mechCoupling',0.11);
bif.setInfluenceFunction(trs.dm.nActuators,2*trs.dm.nActuators-1);
trs.mat.dmIF       = bif.modes;
trs.mat.dmIF_inv = pinv(full(trs.mat.dmIF));
trs.mat.Hdm        = trs.mat.dmIF*trs.mat.dmIF_inv;

%2.5. Tip-tilt
trs.tipTilt.tilt2meter = 1;%!!!!!!
trs.tipTilt.slopes = trs.tipTilt.tilt2meter*aoSys.loopData.tiptilt;
trs.tipTilt.slopes= bsxfun(@minus,trs.tipTilt.slopes,mean(trs.tipTilt.slopes,2));

trs.tipTilt.com = trs.tipTilt.tilt2meter*aoSys.loopData.tiltCom;
trs.tipTilt.com = bsxfun(@minus,trs.tipTilt.com,mean(trs.tipTilt.com,2));
trs.tipTilt.nExp = size(trs.tipTilt.slopes,2);    

%% 3\ Get system matrices
%3.1\ Get DM commands reconstructors from slopes
MC              = aoSys.matrices.DMTTRem*aoSys.wfs2dm.M*aoSys.matrices.SlopeTTRem; %command matrix
trs.mat.R    = trs.dm.volt2meter*MC;
trs.mat.Rtt = trs.tipTilt.tilt2meter;

%3.2\ Get the reconstructed wavefront in OPD and in the actuators space
trs.rec.res    = trs.dm.volt2meter*trs.mat.R*trs.wfs.slopes;
trs.rec.res    = bsxfun(@minus,trs.rec.res,mean(trs.rec.res,2));

%3.3 fill vector to get 21x21 actuators
u = zeros(trs.dm.nActuators^2,trs.wfs.nExp);
u(trs.dm.validActuators,:) = trs.rec.res;
trs.rec.res = u;
u = zeros(trs.dm.nActuators^2,trs.wfs.nExp);
u(trs.dm.validActuators,:) = aoSys.matrices.DMTTRem*trs.dm.com;
trs.dm.com = u;

%% 4\ Get the loop status and model transfer function
%4.1. Latency
if ~isfield(aoSys.parm.loopStatus,'tt')
    aoSys.loopStatus.tt = aoSys.loopStatus.ho;
end
trs.holoop.lat = aoSys.loopStatus.ho.latency;
trs.ttloop.lat = aoSys.loopStatus.tt.latency;

%4.2. Frequency
trs.holoop.freq = 1/aoSys.loopStatus.ho.frameRate/aoSys.tel.samplingTime;
trs.ttloop.freq = 1/aoSys.loopStatus.tt.frameRate/aoSys.tel.samplingTime;

%4.3. RTC controller HO loop: pure integrator
trs.holoop.gain = aoSys.loopStatus.ho.gain;
trs.holoop.tf.num=  [trs.holoop.gain 0 0 0];
trs.holoop.tf.den= [-1 0 0];

trs.ttloop.gain = aoSys.loopStatus.tt.gain;
trs.ttloop.tf.num=  [trs.ttloop.gain 0 0 0];
trs.ttloop.tf.den= [-1 0 0];

%% 5\ Restore the PSF
trs.cam.pixelScale = aoSys.parm.cam.pixelScale;
trs.cam.resolution = size(aoSys.psf,1);
trs.cam.fov = trs.cam.resolution*trs.cam.pixelScale/1e3;

trs.cam.image    =  aoSys.psf ;
trs.cam.wavelength = aoSys.sci.wavelength;
trs.cam.samp          = constants.radian2mas*trs.cam.wavelength/trs.tel.Ddm/2/trs.cam.pixelScale;
trs.src.x                  = 0;
trs.src.y                  = 0;
trs.cam.itime          = aoSys.cam.exposureTime*aoSys.tel.samplingTime;
trs.cam.coadds        = 1;

trs.sky = psfStats( trs.cam.image,trs.tel.pupil,trs.cam.wavelength,...
    trs.cam.samp,trs.cam.pixelScale,trs.src,'flagMoffat',false,'flagGaussian',false);
end