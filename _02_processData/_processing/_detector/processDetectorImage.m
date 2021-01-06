function [trs,flagStatus] = processDetectorImage(trs,varargin)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));
inputs.addParameter('interpBadPix',true,@islogical);
inputs.addParameter('centerOnCog',false,@islogical);
inputs.addParameter('flagGaussian',true,@islogical);
inputs.addParameter('flagMoffat',true,@islogical);
inputs.addParameter('flagBrightestStar',false,@islogical);
inputs.addParameter('ron',0,@isnumeric);
inputs.addParameter('umax',5,@isnumeric);
inputs.addParameter('x0',[],@isnumeric);
inputs.parse(trs,varargin{:});

flagInterpBadPix = inputs.Results.interpBadPix;
centerOnCog      = inputs.Results.centerOnCog;
flagGaussian     = inputs.Results.flagGaussian;
flagMoffat       = inputs.Results.flagMoffat;
flagBrightestStar= inputs.Results.flagBrightestStar;
ron              = inputs.Results.ron;
umax             = inputs.Results.umax;
x0               = inputs.Results.x0;

%% 1\ Read data
path_im = trs.path_imag;
rawImg  = fitsread(path_im);
nCCD    = size(rawImg);
nPSF    = trs.cam.resolution;
if isempty(nPSF)
   nPSF = nCCD/2;
end

%% 2\ Post-processing - step 1: background/flat and bad pixel

%2.1 Crop the calibrated maps if necessary
if all(size(trs.cam.badPixelMap) == nCCD)
    badPixMap = trs.cam.badPixelMap;
elseif all(size(trs.cam.badPixelMap) < nCCD) %cropped mode (not binning)
    badPixMap = puakoTools.crop(trs.cam.badPixelMap,nCCD);   
else % how can we manage this case ?
    badPixMap = 0*rawImg;
end
 
if all(size(trs.cam.flat) == nCCD)
    flatMap = trs.cam.flat;
elseif all(size(trs.cam.flat) < nCCD) %cropped mode (not binning)
    flatMap = puakoTools.crop(trs.cam.flat,nCCD);
else % how can we manage this case ?
    flatMap = ones(nCCD);
end

if all(size(trs.cam.background) == nCCD)
    bgMap = trs.cam.background;
elseif all(size(trs.cam.background) < nCCD) %cropped mode (not binning)
    bgMap = puakoTools.crop(trs.cam.background,nCCD);
else % how can we manage this case ?
     bgMap = 0;
end

if all(size(trs.cam.dark) == nCCD)
    dMap = trs.cam.dark;
elseif all(size(trs.cam.dark) < nCCD) %cropped mode (not binning)
    dMap = puakoTools.crop(trs.cam.dark,nCCD);
else % how can we manage this case ?
    dMap = 0;
end

%2.2. Unbiasing from dark
rawImg = rawImg - dMap;

%2.3. Put bad pixels to zero
badPixMap(rawImg~=rawImg) = 1; % get rid on NAN values
badPixMap(rawImg == 0) = 1;
   
%2.4 Debiasing from the flat field
mskFlat         = flatMap>0.5;
rawImg(mskFlat) = rawImg(mskFlat)./flatMap(mskFlat);
rawImg          = rawImg.*mskFlat;

%2.5. Removing hot pixels
rawImg = removeHotPixelsMap(rawImg,'segmentSize',round(min(size(rawImg))*0.5),'threshold',3);
rawImg = rawImg.*(rawImg==rawImg);

%2.6 Subtract a background
if ~isempty(bgMap) && any(bgMap(:))
    if numel(size(bgMap)) == 3
        bgMap =  median(bgMap,3);
        flagStatus = 'ok';
    elseif  numel(size(bgMap)) == 2
       flagStatus = 'ok';
    else
        uiwait(errordlg('The background map format does not comply the frame format.'));
        flagStatus = 'error';
        return;
    end
else
    % Remove the best bg estimates from the histogram max value
    [val,edges] = histcounts(rawImg,1e4);
    idx = find(val == max(val));
    bgMap = 0.5*(edges(idx) + edges(idx+1));
end
rawImg = rawImg - bgMap;
rawImg(rawImg < -5*std(rawImg)) = 0;
            
%2.7 remove final BG model using a quadratic model
if all(nCCD == 1024)
    res = poly2d_projection(rawImg,2,'msk',badPixMap == 0);
    img = rawImg - res.model;
else
    img = rawImg;
end


%2.8 Interpolating bad pixels from the avalaible neighbors
if any(badPixMap(:)) && flagInterpBadPix    
    supBad = puakoTools.createDeadPixFrame(badPixMap);
    img    = puakoTools.corrDeadPixFrame(supBad, img.*(badPixMap==0) );
end

%2.9 Remove negative pixels
img(img< -5*std(img(:))) = 0;

%% 3\ Detector characteristics and sampling
list                = trs.fitsHdr(:,1);
val                 = trs.fitsHdr(:,2);
trs.cam.wavelength  = cell2mat(val(strcmp(list,'CENWAVE')))*1e-6;
trs.cam.mode        = cell2mat(val(strcmp(list,'CAMNAME')));
switch trs.cam.mode %https://www2.keck.hawaii.edu/inst/nirc2/genspecs.html
    case 'narrow'
        trs.cam.pixelScale = 1e3*0.009942;                                  %pixel scale in mas   
    case 'medium'
        trs.cam.pixelScale = 1e3*0.019829 ;
    case 'wide'
        trs.cam.pixelScale = 1e3*0.039686;
end
trs.cam.fov      = trs.cam.resolution*trs.cam.pixelScale/1e3;        %Field of view in arcsec
trs.cam.samp     = constants.radian2mas*trs.cam.wavelength/trs.tel.Ddm/2/trs.cam.pixelScale;
trs.cam.itime    = cell2mat(val(strcmp(list,'ITIME')));        % integration time
trs.cam.coadds   = cell2mat(val(strcmp(list,'COADDS')));
trs.cam.saturation = 1e4*trs.cam.coadds;

%% 4\ Post-processing - step 2: Stars detection

%1\ Detect stars
[x0_all,y0_all,F0_all] = detectPeaks(img,'border',20,'threshold',3,'iMax',trs.cam.saturation,'flagBrightestStar',flagBrightestStar);

%2\ Count the number of sources and stars in each isoplanatic patches
nSrc = numel(x0_all);
nObj = cellfun(@numel,x0_all);
trs.src = [];
    
%% 5\ Post-processing - step 2: PSF extraction

%1\ Crop the image
trs.cam.image = zeros(trs.cam.resolution,trs.cam.resolution,nSrc);

for iSrc = 1:nSrc
    % unpack stellar parameters
    y0 = y0_all{iSrc};
    x0 = x0_all{iSrc};
    F0 = F0_all{iSrc};
    % select the brightest star
    id = find(F0 == max(F0));
    x0 = x0(id);
    y0 = y0(id);
    % update the structure
    trs.src(iSrc).x    = x0_all{iSrc};
    trs.src(iSrc).y    = y0_all{iSrc};
    trs.src(iSrc).F    = F0_all{iSrc}';
    trs.src(iSrc).nObj = nObj(iSrc);
    % cropping
    if centerOnCog
        % Get the image around the max position with twice resolution
        idx = max(1,round(x0 - nPSF + 1)) : min(nCCD(1),round(x0+nPSF));
        idy = max(1,round(y0 - nPSF + 1)) : min(nCCD(2),round(y0+nPSF));
        img = img(idx,idy);
        % Get the barycenter position
        [~,~,ron]   = puakoTools.getFlux(tmp);
        [xx, yy]    = puakoTools.cog(tmp, ron);
        x0          = round(xx);
        y0          = round(yy);
    end
    idx = max(1,round(x0 - nPSF/2 + 1)) : min(nCCD(2),round(x0+nPSF/2));
    idy = max(1,round(y0 - nPSF/2 + 1)) : min(nCCD(1),round(y0+nPSF/2));
    img_i = img(idy,idx);
    
    %remove bad pixels
    
    %Zero-pad if the sources are too close from edges
    dx = 0;
    dy = 0;
    if size(img_i,1)<trs.cam.resolution
        if idy(1) == 1
            img_i = padarray(img_i,[trs.cam.resolution-size(img_i,1),0],'pre');
            dy =  -(trs.cam.resolution-size(img_i,1));
        else
            img_i = padarray(img_i,[trs.cam.resolution-size(img_i,1),0],'post');
            dy =  (trs.cam.resolution-size(img_i,1));
        end
    end
    if size(img_i,2)<trs.cam.resolution
        if idx(1) == 1
            img_i = padarray(img_i,[0,trs.cam.resolution-size(img_i,2)],'pre');
            dx =  -(trs.cam.resolution-size(img_i,2));
        else
            img_i = padarray(img_i,[0,trs.cam.resolution-size(img_i,2)],'post');
            dx =  (trs.cam.resolution-size(img_i,2));
        end
    end
    % allocation
    trs.cam.image(:,:,iSrc) = img_i;
end

%2\ Update the sources positions

if nSrc > 1
    if strcmp(trs.aoMode,'NGS')
        % In NGS mode, the guide star is visible in the image
        % We assume the NGS is the brightest star in the fov
        Fmax = max([trs.src.F]);
        idG  = zeros(1,nSrc);
        for k = 1:nSrc
            idG(k) = any(trs.src(k).F == Fmax);
        end
        idG = find(idG);
        trs.ngs.x = (trs.src(idG).x - nCCD(2)/2 - 1)*trs.cam.pixelScale/1e3;
        trs.ngs.y = (trs.src(idG).y - nCCD(1)/2 - 1)*trs.cam.pixelScale/1e3;
        for iSrc = 1:nSrc
            trs.src(iSrc).x = (x0_all{iSrc} - nCCD(2)/2 - 1)*trs.cam.pixelScale/1e3 - trs.ngs.x;
            trs.src(iSrc).y = (y0_all{iSrc} - nCCD(1)/2 - 1)*trs.cam.pixelScale/1e3 - trs.ngs.y;
        end
    else
        %In LGS mode, we do not see the laser due to the dichroic and the NIR filter
        % We assume the LGS is pointing the center of the image
        trs.lgs.x = 0;
        trs.lgs.y = 0;
        for iSrc = 1:nSrc
            trs.src(iSrc).x = (trs.src(idG).x - nCCD(2)/2 - 1)*trs.cam.pixelScale/1e3;
            trs.src(iSrc).y = (trs.src(idG).y - nCCD(1)/2 - 1)*trs.cam.pixelScale/1e3;
        end
    end
else
    trs.ngs.x = 0;trs.ngs.y = 0;
    trs.lgs.x = 0;trs.lgs.y = 0;
    trs.src.x = 0;trs.src.y = 0;
end

%% 6\ Post-processing - step 3: PSF metrics estimation
trs.sky = psfStats.empty(nSrc,0);
for iSrc = 1:nSrc
    trs.sky(iSrc) = psfStats(trs.cam.image(:,:,iSrc),trs.tel.pupil,trs.cam.wavelength,...
        trs.cam.samp,trs.cam.pixelScale,trs.src(iSrc),'ron',ron,...
        'flagMoffat',flagMoffat,'flagGaussian',flagGaussian,'umax',umax,'x0',x0);
end

%% 7\ Differential field dependent aberrations
path_stat = [trs.path_calibration,trs.cam.name,'/STATIC/'];
folder_stat = dir(path_stat);
idxstat = find(contains(upper({folder_stat.name}),'MAP') & contains(upper({folder_stat.name}),'.FITS'));
if ~isempty(idxstat)
    trs.cam.field_static_map    = fitsread([path_stat,folder_stat(idxstat).name]);
    trs.cam.x_stat              = trs.cam.pixelScale*(fitsread([path_stat,'x_vals.fits']) - 511)/1e3;
    trs.cam.y_stat              = trs.cam.pixelScale*(fitsread([path_stat,'y_vals.fits']) - 511)/1e3;
    trs.cam.diff_field_stat     = differentialStaticAberrations(trs.src,trs.cam);
end
