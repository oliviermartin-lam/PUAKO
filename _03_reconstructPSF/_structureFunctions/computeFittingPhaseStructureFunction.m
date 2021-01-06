%{
------------HEADER-----------------
Objective          ::  Compute bi-dimensional phase structure function of the
aliasing phase

INPUT VARS
 psfr          :: The psfReconstruction top-level class
 
OUTPUT VARS
 sf_2D             :: bi-dimensional phase structure function map (Toeplitz) of the fitting error (psfr.nOtf x psfr.nOtf)
psd                  :: Normalized (r0 = 1 m) Power spectrum density of the fitting error (psfr.nOtf x psfr.nOtf)
Created by         :: O. Beltramo-Martin - ONERA/LAM
Creation date      :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}
function obj = computeFittingPhaseStructureFunction(obj,varargin)
inputs = inputParser;
inputs.addRequired('obj',@(x) isa(x,'psfReconstruction'));
inputs.addParameter('aoPattern','influ',@ischar)
inputs.parse(obj,varargin{:});
aoPattern = inputs.Results.aoPattern;

%1\ Parsing inputs
r0 = obj.trs.res.seeing.r0*(obj.trs.cam.wavelength/0.5e-6)^1.2;
L0 = obj.trs.res.seeing.L0;
nActu = obj.trs.dm.nActuators;
nT = obj.otf.nTimes;
d = obj.trs.dm.pitch;
%2\ Define the frequency space
nK   = obj.otf.nOtf;
[kx,ky] = freqspace(nK,'meshgrid');
kc  = 1/(2*d);
kx = kx*kc*nT;
ky = ky*kc*nT;
k   = hypot(kx,ky);
dk  = 2*kc*nT/nK; % Pixel scale

%3\ Define the atmospheric phase PSD
cst = r0^(-5/3)*(24*gamma(6/5)/5)^(5/6)*(gamma(11/6)^2/(2*pi^(11/3)));
psd = cst*(k.^2 + 1/L0.^2).^(-11/6);

%4\ Filtering the AO-controlled area

if strcmp(aoPattern,'circle')
    msk = k > kc;
elseif strcmp(aoPattern,'square')
    msk = ~(abs(kx) <= kc & abs(ky) <= kc);
else
    msk = 1-obj.trs.dm.modes.getTransferFunction(21,nK);
    msk(msk<1e-2) = 0;
end

%5\ Update the PSFR class
obj.cov.psdFit= psd.*msk;
obj.sf.Dfit = puakoTools.cov2sf(puakoTools.psd2cov(obj.cov.psdFit,dk));
end


