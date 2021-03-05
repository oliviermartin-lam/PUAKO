parm = [];
if ~exist('atmProfile','var')
    atmProfile = 'median';
end
atmProfile = 'sparse';
%% ATM
if ~exist('zenithAngle','var')
    parm.atm.zenithAngle        = 0;         % rad
else
    parm.atm.zenithAngle        = zenithAngle;         % rad
end
parm.atm.zenithAngle            = pi/6;         % rad
airmass                         = 1/cos(parm.atm.zenithAngle);

parm.atm.photometry             = photometry.V0;% photometric band
parm.atm.windSpeed              = 10;           % wind speed vector in m/s
parm.atm.windDirection          = 0;            % wind direction vector in rad
parm.atm.r0                     = 0.16*cos(parm.atm.zenithAngle )^(3/5); % coherence lenght in meters at 0.5microns
parm.atm.L0                     = 50;           % Outer scale in meters
parm.atm.altitude               = [0 0.5 1 2 4 8 16].*1e3./cos(parm.atm.zenithAngle ); % altitude vector in meters
parm.atm.windSpeed              = [6.8 6.9 7.1 7.5 10.0 26.9 18.5];
parm.atm.windDirection          = [0 pi/2 pi/4 3*pi/2 8*pi/3 pi/8 pi];    

% fractional Cn2 profile (sum=1)
switch atmProfile
    case 'median'
        parm.atm.fractionalR0  = [0.517 0.119 0.063 0.061 0.105 0.081 0.054]; %median
    case 'challenging'
        parm.atm.fractionalR0  = [0.464 0.109 0.060 0.062 0.135 0.113 0.057]; %challenging
    case 'sparse'
        parm.atm.fractionalR0  = [0.7 0.3]; %challenging
        parm.atm.altitude      = [0.0 1e4];
        parm.atm.windSpeed     = [6.8,15.0];
        parm.atm.windDirection = [0 pi/4];
end




% %%
% parm = [];
% 
% %% ATM
% parm.atm.zenithAngle            = pi/6;                         % rad
% airmass                         = 1/cos(parm.atm.zenithAngle);
% parm.atm.photometry             = photometry.V0;                % photometric band
% parm.atm.windSpeed              = [10,25,15];                   % wind speed vector in m/s
% parm.atm.windDirection          = [0,pi/3,pi/2];                % wind direction vector in rad
% parm.atm.r0                     = 0.12/airmass^(3/5);           % coherence lenght in meters at 0.5microns
% parm.atm.L0                     = 25;                           % Outer scale in meters
% parm.atm.altitude               = [0,5,10]*1e3*airmass;         % altitude vector in meters
% parm.atm.fractionalR0           = [.5,.3,.2];                   % fractional Cn2 profile (sum=1)

%% SCIENCE                                              
parm.sci.magnitude              = 0;                            % stellar magnitude
parm.sci.photometry             = photometry.K;                 % photometric band
parm.sci.height                 = Inf;                          % height (default Inf)
parm.sci.x                      = 0;                            % Sources x cartesian position in arcsec
parm.sci.y                      = 0;                            % Sources x cartesian position in arcsec

%% TEL
parm.tel.D                      = 11.25;                        % telescope primary mirror diameter
parm.tel.obstructionRatio       = 0.2356;                       % central obscuration ratio
parm.tel.fieldOfViewInArcsec    = 10;                           % fieldOfViewInArcsec
parm.tel.resolution             = 200;                          % resolution 
parm.tel.samplingTime           = 1e-3;                         % sampling time in sec
parm.tel.pupilShape             = [];                           % User-defined pupil-shape

%% SEGMENTED PUPIL
parm.seg.usePupilClass          = true;
parm.seg.nSides                 = 6;                            % Number of segment sides
parm.seg.radius                 = 0.9;                          % Segment radius
parm.seg.resolution             = parm.tel.resolution;         % Segment resolution in pixel
xSeg                            = parm.seg.radius*[1.5 0 -1.5 -1.5 0 1.5 3 1.5 0 -1.5 -3 -3 -3 -1.5 0  1.5 3 3 4.5 3 1.5 0 -1.5 -3 -4.5 -4.5 -4.5 -4.5 -3 -1.5 0 1.5 3 4.5 4.5 4.5];
ySeg                            = sqrt(3)/2*parm.seg.radius*[1 2 1 -1 -2 -1  2  3  4   3   2  0 -2  -3 -4   -3 -2 0 3 4 5 6 5 4 3 1 -1 -3 -4 -5 -6 -5 -4 -3 -1 1];
parm.seg.segmentCoordinates     = [xSeg' ySeg'];                % Segment coordinates
parm.seg.pupilUnit              = 'm';                          % Units for coordinates: pixel | meter
parm.seg.pupilSpider.n          = 3;                            % Number of spiders
parm.seg.pupilSpider.angle      = [0 pi/3 2*pi/3];              % Spiders angular positions
parm.seg.pupilSpider.width      = 0.1;                        % Spider width in pupilUnit
parm.seg.pupilSpider.usrDefined = 0;
parm.seg.pupilSpider.usrMask    = [];

diffPiston                      = 0*100e-9;                               % Statistical diff piston in m;
rampPiston                      = 0*80e-9;                                % Piston ramp max value in m                  
piston                          = diffPiston*randn(1,length(parm.seg.segmentCoordinates));
piston                          = piston-mean(piston);
ramp                            = keckTools.keck_pistonRamp(rampPiston);
parm.seg.pupilCoeffPhaseModes   = piston' + ramp;               % Vector of piston errors
parm.seg.pupilReflexivity       = [];                           % Vector of transmission disparities


%% LGS SOURCES
parm.lGs.magnitude              = 11;                           % magnitude                         
parm.lGs.photometry             = photometry.Na;                % photometric band                                       
parm.lGs.height                 = 90e3*airmass;                 % laser beacon mean height    
parm.lGs.nPhoton                = 1000*20^2/10^2*1e3;           % number of photons/m^2/s -> 1000 ph/frames/subap;
parm.lGs.naProfile              = ones(1,length(parm.lGs.height));% vector with normalised profile strenghts
parm.lGs.spotFwhm               = [0];                           % laser spot FWHM                    
parm.lGs.x                      = [0];                      % LGSs x cartesian position in arcsec
parm.lGs.y                      = [0];                     % LGSs y cartesian position in arcsec

%% NGS SOURCES
parm.nGs.magnitude              = 11;                           % magnitude                         
parm.nGs.photometry             = photometry.R;                 % photometric band                    
parm.nGs.x                      = 0;                            % NGSs x cartesian position in arcsec
parm.nGs.y                      = 0;                            % NGSs y cartesian position in arcsec                   

%% WAVE-FRONT SENSOR
parm.wfs.nLenslet               = 20;                           % number of lenslets
parm.wfs.nPx                    = 10;                           % number of pixels per lenslet
parm.wfs.pixelScale             = parm.wfs.nLenslet*constants.radian2mas*parm.nGs.photometry.wavelength/parm.tel.D/2;          % Pixel scale in mas          
parm.wfs.minLightRatio          = 0.5;                          % [0-1] partial illumination ratio to consider valid a sub-aperture
parm.wfs.exposureTime           = 1e-3;                         % WFS exposure time
parm.wfs.frameRate              = 1;                            % slaved to tel.samplingTime
parm.wfs.ron                    = 0;                            % ron in e-
parm.wfs.photonNoise            = false;                        % true/false statement
parm.wfs.pixelThreshold         = parm.wfs.ron/2;               % Pixel threshold in the cog algorithm
parm.wfs.throughput             = .6;                           % [0-1] throughput ratio 

%% TT WFS
parm.lowfs.nLenslet             = 1;                            % number of lenslets
parm.lowfs.nPx                  = 100;                          % number of pixels per lenslet
parm.lowfs.pixelScale           = constants.radian2mas*parm.nGs.photometry.wavelength/parm.tel.D;          % pixel scale
parm.lowfs.minLightRatio        = 0;                            % [0-1] partial illumination ratio to consider valid a sub-aperture
parm.lowfs.exposureTime         = 1e-3;                         % WFS sampling time
parm.lowfs.frameRate            = 1;                            % slaved to tel.samplingTime
parm.lowfs.ron                  = 0;                          % ron in e-
parm.lowfs.photonNoise          = false;                         % true/false statement
parm.lowfs.pixelThreshold       = 0;                            % Pixel threshold in the cog algorithm
parm.lowfs.throughput           = 0.6;                          % [0-1] throughput ratio 

%% DEFORMABLE MIRROR
parm.dm.pitch                   = 0.5625;                          % dm pitch
parm.dm.nActuators              = 21;                           % Number of actuators par row
parm.dm.crossCoupling           = 0.2;                          % actuator cross-coupling
parm.dm.influenceType           = 'gaussian';                   % type of influence function: gaussian | monotonic | overshoot | other

%% IMAGING CAMERA
parm.cam.pixelScale             = 9.94;                         % Pixel scale in mas                         
parm.cam.resolution             = 200;                          % camera resolution                               
parm.cam.clockRate              = 1;                            % clock-rate (default 1)                   
parm.cam.exposureTime           = 5000;                         % exposure time in time-steps x clockRate                   
parm.cam.startDelay             = 10;                           % transient period for loop to converge in time-steps x clockRate                   
parm.cam.ron                    = 0;                            % ron in e-                  
parm.cam.photonNoise            = false;                        % true/false photon noise   

%% LOOPS SETTINGS
parm.loopStatus.ho.nModeFiltered= 5;                            % # modes to truncate during SVD
parm.loopStatus.ho.controlBasis = [];                           % controlBasis: classic (default) | modal
parm.loopStatus.ho.gain         = 0.5;                          % integrator control loop gain
parm.loopStatus.ho.latency      = 1e-3;                         % loop latency (default tel.samplingTime)
parm.loopStatus.ho.tomography   = [];                           % tomographic algorithm
parm.loopStatus.ho.control      = 'integrator';                 % controller
parm.loopStatus.ho.gainPolc     = 1;                            % POLC control loop gain

parm.loopStatus.lo.gain         = 0.5;                          % integrator control loop gain       
parm.loopStatus.lo.latency      = 1e-3;                         % loop latency (default tel.samplingTime)
parm.loopStatus.lo.tomography   = [];                           % tomographic algorithm: mean | mmse
parm.loopStatus.lo.control      = 'integrator';                 % controller : integrator | polc
parm.loopStatus.lo.gainPolc     = 1;                            % POLC control loop gain
%% NCPA
parm.ncpa                       = [];                           % NCPA map in meter




