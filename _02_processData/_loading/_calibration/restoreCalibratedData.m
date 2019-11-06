function trs = restoreCalibratedData(trs)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));

date = trs.date;
path_calib = trs. path_calibration;

%% 1\ Get static aberrations calibration

if strcmp(date,'20130801')
    trs.tel.static_map = fitsread([path_calib,'_phaseDiversity/phasemaps_2013.fits'])*1.6455e3/2/pi;
    trs.tel.static_map = squeeze(trs.tel.static_map(:,:,2));
else
    pc       = [path_calib,'_phaseDiversity/phase_diversity_results_',trs.date,'_average_PD1.sav'];
    cal      = restore_idl('filename',pc);
    ncpaMap1 = double(cal.PHASE);
    ncpaWvl  = double(cal.LAMBDA)*1e-6;
    pc(end-4)= '2';
    cal      = restore_idl('filename',pc);
    ncpaMap  = double(cal.PHASE)*0.5 + 0.5*ncpaMap1;
    trs.tel.static_map = tools.rotateIm(ncpaMap,90)*ncpaWvl*1e9/2/pi;
end

%% 2\ Get the pupil shape
trs.tel.resolution = size(trs.tel.static_map,1); 
trs.tel.pixelscale = trs.tel.Ddm/trs.tel.resolution;  
trs.cam.pupilTracking = upper(cell2mat(trs.fitsHdr.value(strcmp(trs.fitsHdr.field,'PMSNAME'))));
trs.tel.elevation = cell2mat(trs.fitsHdr.value(strcmp(trs.fitsHdr.field,'EL')));
trs.tel.azimuth = cell2mat(trs.fitsHdr.value(strcmp(trs.fitsHdr.field,'AZ')));
trs.tel.zenith_angle = 90 - trs.tel.elevation;
trs.tel.airmass = 1/cos(trs.tel.zenith_angle*pi/180);
trs.tel.pupilAngle  = cell2mat(trs.fitsHdr.value(strcmp(trs.fitsHdr.field,'ROTPPOSN'))) - trs.tel.elevation; %See WFS_RotationAngle.pdf - C. Neyman

switch trs.cam.pupilTracking        
    case 'LARGEHEX'
        tmp = restore_idl([path_calib,'_pupil/LargeHex.sav']);       
        [X,Y] = getGridCoordinates(trs.tel.resolution,trs.tel.resolution,tmp.B);
        outerHex = hexagonalSegment(X,Y,2*tmp.B,trs.tel.pupilAngle);
        innerHex = hexagonalSegment(X,Y,2*tmp.A,trs.tel.pupilAngle);
        trs.tel.pupil = outerHex - innerHex;
    case 'INCIRCLE'
        tmp = restore_idl([path_calib,'_pupil/KeckSegmentedPupilArchitecture.sav']);
        dM1 = trs.tel.pupil.Dcircle;
        dM2 = double(tmp.MIR.DINT);
        [~,~,R] = getGridCoordinates(trs.tel.resolution,trs.tel.resolution,dM1/0.968/2);
        trs.tel.pupil = (R<=dM1/2).*(R>dM2/2);        
    otherwise
        if ~isempty(trs.tel.static_map) && any(trs.tel.static_map(:))
            trs.tel.pupil = double(logical(trs.tel.static_map));
        else
            trs.tel.pupil = tools.interpolate(fitsread([path_calib,'_pupil/keckPupil.fits']),trs.tel.resolution);
        end
end


%% 2\ Get DM calibration

%1\ Valid apertures/actuators
trs.dm.validActuators     = logical(load([path_calib,'_dm/keckValidActuatorMap.txt']));
trs.wfs.validSubaperture= logical(load([path_calib,'_dm/keckValidSubapMap.txt']));

%2\ Influence DM functions
bif      = xineticsInfluenceFunctionModel(trs.dm.pitch,'mechCoupling',0.11);
bif.setInfluenceFunction(trs.dm.nActuators,2*trs.dm.nActuators-1);

trs.mat.dmIF       = bif.modes;
trs.mat.dmIF_inv = pinv(full(trs.mat.dmIF));
trs.mat.Hdm        = trs.mat.dmIF*trs.mat.dmIF_inv;

%% 3\ Detector data

if strcmp(date,'20130801')
    badPixMap          = fitsread([path_calib,'_badPixels/supermask2013_modified.fits']);
    badPixMap(169,179) = 1;
    badPixMap(39,707)  = 1;
    flat               = fitsread([path_calib,'_flat/flat_kp2017.fits']);
else
    badPixMap = fitsread([path_calib,'_badPixels/supermask2017.fits']);
    flat      = fitsread([path_calib,'_flat/flat_kp2017.fits']);
    badPixMap(289,136) = 1;
    badPixMap(290,136) = 1;
    badPixMap(289,137) = 1;
    badPixMap(290,137) = 1;
    badPixMap(694,860) = 1;
    badPixMap(578,924) = 1;
    badPixMap(720,859) = 1;
    badPixMap(411,762) = 1;
    badPixMap(449,747) = 1;
end

trs.cam.badPixelMap = badPixMap;                    
trs.cam.background = 0;                                    
trs.cam.flat = flat;            

