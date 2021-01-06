function trs = restoreCalibratedData(trs,getImageOnly)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));
inputs.addRequired('getImageOnly',@ischar);

date = trs.date;
path_calib = trs. path_calibration;
trs.cam.name = cell2mat(trs.fitsHdr(strcmp(trs.fitsHdr(:,1),'CURRINST'),2));

path_calib_cam = [path_calib,trs.cam.name,'/'];
path_calib_ao = [path_calib,'AOSYSTEM/'];
path_calib_tel = [path_calib,'TELESCOPE/'];

%% 1\ Get static aberrations calibration
%1.1 NCPA
if isempty(trs.path_ncpa)
    if strcmp(date(1:4),'2013')
        trs.path_ncpa = [path_calib_ao,'STATIC/phasemaps_2013.fits'];
        trs.tel.static_map = fitsread(trs.path_ncpa)*1.6455e3/2/pi;
        if strcmp(date(5:6),'02')
            trs.tel.static_map = squeeze(trs.tel.static_map(:,:,1));
        elseif strcmp(date(5:6),'08')
            trs.tel.static_map = squeeze(trs.tel.static_map(:,:,2));
        elseif strcmp(date(5:6),'09')
            trs.tel.static_map = squeeze(trs.tel.static_map(:,:,3));
        end
        trs.cam.x_ncpa      = 0; %???
        trs.cam.y_ncpa      = 0; % ???
        cal                 = [];
        cal.lambda          = 1.6455;
        WFSTHETA = 0;
    elseif strcmp(date(1:4),'2017')
        trs.path_ncpa       = [path_calib_ao,'STATIC/phase_diversity_results_',trs.date,'_average_PD2.sav'];
        cal                 = restore_idl('filename',trs.path_ncpa);
        ncpaWvl             = double(cal.LAMBDA)*1e-6;
        trs.tel.static_map  = double(cal.PHASE)*ncpaWvl*1e9/2/pi;
        trs.cam.x_ncpa      = 0; %???
        trs.cam.y_ncpa      = 0; % ???
        WFSTHETA            = cal.WFSTHETA;
    end
else
    if contains(trs.path_ncpa,'.sav')
        cal                 = restore_idl('filename',trs.path_ncpa);
        ncpaWvl             = double(cal.LAMBDA)*1e-6;
        trs.tel.static_map  = double(cal.PHASE)*ncpaWvl*1e9/2/pi;
        trs.cam.x_ncpa      = 0; %???
        trs.cam.y_ncpa      = 0; % ???
        WFSTHETA            = cal.WFSTHETA;
    else
        trs.tel.static_map  = mean(fitsread(trs.path_ncpa),3)*1.6455e3/2/pi; 
        WFSTHETA            = 0;
    end
end

%1.2 Static fitting error
folder_stat = dir([path_calib_ao,'STATIC/']);
idx         = contains(upper({folder_stat.name}),'FITTING') & contains(upper({folder_stat.name}),'.FITS');
if any(idx)
    stat_file               = {folder_stat(idx).name};
    idx                     = find(idx);
    nStat                   = numel(stat_file);
    trs.tel.static_fitting  = 0;
    for k=1:nStat
        trs.tel.static_fitting  = trs.tel.static_fitting + 410*fitsread([path_calib_ao,'STATIC/',folder_stat(idx(k)).name])/nStat; %in nm
    end
    trs.tel.static_fitting      = puakoTools.rotateIm(trs.tel.static_fitting,90);
end

%% 2\ Get the pupil shape

trs.tel.resolution      = 512; %the final pupil resolution will be set up during the reconstruction process
trs.tel.pixelscale      = trs.tel.Ddm/trs.tel.resolution;
trs.cam.pupilTracking   = upper(cell2mat(trs.fitsHdr(strcmp(trs.fitsHdr(:,1),'PMSNAME'),2)));
trs.tel.elevation       = cell2mat(trs.fitsHdr(strcmp(trs.fitsHdr(:,1),'EL'),2));
trs.tel.azimuth         = cell2mat(trs.fitsHdr(strcmp(trs.fitsHdr(:,1),'AZ'),2));
trs.tel.zenith_angle    = 90 - trs.tel.elevation;
trs.tel.airmass         = 1/cos(trs.tel.zenith_angle*pi/180);

% this part must be verified !

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
theta = (cell2mat(trs.fitsHdr(strcmp(trs.fitsHdr(:,1),'ROTPPOSN'),2)) - trs.tel.elevation -  cell2mat(trs.fitsHdr(strcmp(trs.fitsHdr(:,1),'INSTANGL'),2)));
if theta < 0
    trs.tel.pupilAngle  = 90 - (cell2mat(trs.fitsHdr(strcmp(trs.fitsHdr(:,1),'ROTPPOSN'),2)) - trs.tel.elevation -  cell2mat(trs.fitsHdr(strcmp(trs.fitsHdr(:,1),'INSTANGL'),2))); %See WFS_RotationAngle.pdf - C. Neyman
else
    trs.tel.pupilAngle  =  -theta ;
end
trs.wfs.pupilAngle    = -90 + WFSTHETA -trs.tel.pupilAngle ;
%trs.wfs.pupilAngle = WFSTHETA -trs.tel.pupilAngle ;
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

switch trs.cam.pupilTracking
    case 'LARGEHEX'
        tmp             = restore_idl([path_calib_tel,'LargeHex.sav']);
        res             = getGridCoordinates(trs.tel.resolution,trs.tel.resolution,tmp.B);
        osamp           = round(1500/trs.tel.resolution);
        outerHex        = hexagonalSegment(res.x2D,res.y2D,2*tmp.B,trs.tel.pupilAngle,'overSamp',osamp);
        innerHex        = hexagonalSegment(res.x2D,res.y2D,2*tmp.A,trs.tel.pupilAngle,'overSamp',osamp);
        trs.tel.pupil   = double(logical(outerHex - innerHex));
    case 'INCIRCLE'
        tmp             = restore_idl([path_calib_tel,'KeckSegmentedPupilArchitecture.sav']);
        dM1             = trs.tel.Dcircle;
        dM2             = double(tmp.MIR.DINT);
        res             = getGridCoordinates(trs.tel.resolution,trs.tel.resolution,dM1/0.968/2);
        trs.tel.pupil   = (res.r2D<=dM1/2).*(res.r2D>dM2/2);
    otherwise
        if ~isempty(trs.tel.static_map) && any(trs.tel.static_map(:))
            trs.tel.pupil = double(logical(trs.tel.static_map));
        else
            trs.tel.pupil = puakoTools.interpolate(fitsread([path_calib_tel,'keckPupil.fits']),trs.tel.resolution);
        end
end
trs.tel.static_map  = puakoTools.rotateIm(trs.tel.static_map,trs.wfs.pupilAngle);
trs.dm.pupilMask    = logical(puakoTools.interpolate(trs.tel.pupil,trs.dm.nActuators));

    %% 2\ Get DM calibration
if ~getImageOnly
   
    %1\ Valid apertures/actuators
    trs.dm.validActuators       = logical(load([path_calib_ao,'keckValidActuatorMap.txt']));
    trs.wfs.validSubaperture    = logical(load([path_calib_ao,'keckValidSubapMap.txt']));
    
    %2\ Influence DM functions
    bif      = xineticsInfluenceFunctionModel(trs.dm.pitch,'mechCoupling',0.11);
    bif.setInfluenceFunction(trs.dm.nActuators,2*trs.dm.nActuators-1);
    % rotation between IDL and Matlab convention
    % IDL: bottom-left corner to top right, left to right at each row
    %Matlab: top-left corner to bottom right, top to bottom each row
    % 90deg rotation counterclockwise needed
    % Actually not, otherwise the reconstructed phase variance map does not
    % comply with the WFS intensity map (more variance when less light)
    trs.mat.dmIF        = bif.modes;%puakoTools.idlToMatlabIndex(bif.modes,trs.dm.nActuators,true(trs.dm.nActuators),'modes');
    trs.mat.dmIF_inv    = pinv(full(trs.mat.dmIF));
    trs.mat.Hdm         = trs.mat.dmIF*trs.mat.dmIF_inv;
end
%% 3\ Detector data

if isfile(trs.path_imag)
    trs.cam.rawFrame = fitsread(trs.path_imag);
end

%3.1 DARK 
tmp = split(replace(trs.path_imag,'//','/'),'/');
path_img = [strjoin(tmp(1:end-1),'/'),'/'];
folder_img = dir(path_img);
    
% Read the night folder
if isempty(trs.path_dark)    
    idx = find(contains(upper({folder_img.name}),'DARK'));
    is_dark = any(idx);
    
    if is_dark
        trs.path_dark = [path_img,folder_img(idx).name,'/'];
    else
        % Read the calibration folder
        folder_calib = dir(path_calib);
        idx = find(contains(upper({folder_calib.name}),'DARK'));
        trs.path_dark = [path_img,folder_calib(idx).name,'/'];
    end
end
trs.cam.dark = getCalibratedDectorData(trs.path_dark,trs.fitsHdr,{'ITIME','COADDS'});


%3.2 flat field
idx = find(contains(upper({folder_img.name}),'FLAT'));
if ~any(idx)
    path_flat = [path_calib_cam,'FLAT/'];
    folder_flat = dir(path_flat);
    idxflat = find(contains(upper({folder_flat.name}),'FLAT') & contains(upper({folder_flat.name}),'.FITS'));
    if ~isempty(idxflat)
        nF = numel(idxflat);
        trs.cam.flat = 0;
        for k=1:nF
            trs.cam.flat = trs.cam.flat + fitsread([path_flat,folder_flat(idxflat(k)).name])/nF;
        end
    end
else
    trs.cam.flat = 1;
end

%3.3. Background
idx = find(contains(upper({folder_img.name}),'BACK'));
if ~any(idx)
    path_sky = [path_calib_cam,'BACK/'];
    folder_sky = dir(path_sky);
    idxsky = find(contains(upper({folder_sky.name}),'BACK') & contains(upper({folder_sky.name}),'.FITS'));
elseif ~isempty(trs.path_sky)
    folder_sky = dir(trs.path_sky);
    idxsky = find(contains(upper({folder_sky.name}),'BACK') & contains(upper({folder_sky.name}),'.FITS'));
end
if ~isempty(idxsky)
    nF = numel(idxsky);
    trs.cam.background = 0;
    for k=1:nF
        trs.cam.background = trs.cam.background + fitsread([path_sky,folder_sky(idxsky(k)).name])/nF;
    end
end


%3.4 Bad pixels
path_pix = [path_calib_cam,'BADPIXELS/'];
folder_pix = dir(path_pix);

idxpix = find(contains(upper({folder_pix.name}),'MASK') & contains(upper({folder_pix.name}),'.FITS') & contains(upper({folder_pix.name}),trs.date(1:4)));
if ~isempty(idxpix)
    nMap = numel(idxpix);
    trs.cam.badPixelMap = false;
    for k=1:nMap
        trs.cam.badPixelMap = trs.cam.badPixelMap | fitsread([path_pix,folder_pix(idxpix(k)).name]);
    end    
end
trs.cam.badPixelMap(302,749) = 1;
trs.cam.badPixelMap(302,750) = 1;



