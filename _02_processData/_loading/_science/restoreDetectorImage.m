function [trs,flagStatus] = restoreDetectorImage(trs,varargin)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));
inputs.addParameter('interpBadPix',true,@islogical);
inputs.parse(trs,varargin{:});

flagInterpBadPix = inputs.Results.interpBadPix;

%% 1\ Read data
path_im = trs.path_imag;
rawImg  = fitsread(path_im);
nCCD = size(rawImg,1);
nPSF  = trs.cam.resolution;
if isempty(nPSF)
   nPSF = nCCD/2;
end

trs.cam.name = cell2mat(trs.fitsHdr.value(strcmp(trs.fitsHdr.value,'CURRINST')));

%% 2\ Post-processing - step 1: background/flat and bad pixel

%3.1 Scale the badpixel and flat maps if necessary
if size(trs.cam.badPixelMap,1) == nCCD
    badPixMap = trs.cam.badPixelMap;
    flatMap = trs.cam.flat;
    bgMap = trs.cam.background;    
else %cropped mode (not binning)
    badPixMap = tools.crop(trs.cam.badPixelMap,nCCD);
    flatMap = tools.crop(trs.cam.flat,nCCD);
    bgMap = tools.crop(trs.cam.background,nCCD);
end
          
%3.2 Debiasing from the background
if ~isempty(bgMap) && any(bgMap(:))
    if numel(size(bgMap)) == 3
        bgMap =  median(bgMap,3);
        flagStatus = 'ok';
    elseif  numel(size(bgMap)) == 2
       flagStatus = 'ok';
    else
        uiwait(errordlg('The specified path does not refer to any existing folder or file.'));
        flagStatus = 'error';
        return;
    end
else
    bgMap = 0;
end
rawImg = rawImg - bgMap;

%3.3 Debiasing from the flat field
rawImg = rawImg./flatMap;
rawImg(rawImg~=rawImg) = 0; % get rid on NAN values

%3.4 Put bad pixels to zeros           
rawImg = ~badPixMap.*rawImg;

%3.5 Subtract the median value
img = rawImg - median(rawImg(~badPixMap));
            
%3.6 Interpolating bad pixels from the avalaible neighbors
if any(badPixMap(:)) && flagInterpBadPix
    supBad = tools.createDeadPixFrame(badPixMap);
    img    = tools.corrDeadPixFrame(supBad, rawImg );
end
            
%% 3\ Post-processing - step 2: PSF cropping

%3.1 Get rid of image edges
border = 20;
img2 = img(border+1:end-border,border+1:end-border);

%3.2 Perform a median filtering to get rid of remaining hot pixels
n    = 1;
for i=1:n
    img2(:) = medfilt1(img2(:));
end
%3.3 Localize the PSF peak position on the filtered image
[x0,y0] = find(img2 == max(img2(:)));
x0 = x0 + border;
y0 = y0 + border;
            
%3.4. Localize the PSF peak on the real image 
nLoc = length(x0);
pVal = zeros(1,nLoc);
for i=1:nLoc
    pVal(i) = img(x0(i),y0(i));
end
x0 = x0(find(pVal == max(pVal)));
y0 = y0(find(pVal == max(pVal)));
            
%3.5 Crop PSF to the wished field of view
if ~isempty(nPSF) && any(nPSF)
    if trs.cam.resolution < 2*min([nCCD-x0 nCCD-y0])
        idx = round(x0 - nPSF/2 + 1):round(x0+nPSF/2);
        idy = round(y0 - nPSF/2 + 1):round(y0+nPSF/2);
        img = img(idx,idy);
    end
end
               

%% 4\ Fill in the trs structure
list      = trs.fitsHdr.field;
val       = trs.fitsHdr.value;
trs.cam.wavelength = cell2mat(val(strcmp(list,'CENWAVE')))*1e-6;
trs.cam.samp           = 1.5*(trs.cam.wavelength/1.6455e-6);
trs.cam.frame          = img;
trs.src.x                   = (x0-nCCD/2)*trs.cam.pixelScale/1e3;
trs.src.y                   = (y0-nCCD/2)*trs.cam.pixelScale/1e3;

if strcmp(trs.aoMode,'NGS')
    trs.ngs.x = trs.src.x;
    trs.ngs.y = trs.src.y;
else
    trs.lgs.x = trs.src.x;
    trs.lgs.y = trs.src.y;
end

