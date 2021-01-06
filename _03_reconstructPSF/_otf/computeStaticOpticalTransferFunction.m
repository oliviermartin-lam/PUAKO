%{
------------HEADER-----------------
Objective       ::  Compute bi-dimensional phase structure function of the aliasing phase

INPUT VARS
psfr            :: The psfReconstruction top-level class
addStatModes    ::
flagStaticMap   ::

OUTPUT VARS
otfDL           :: Optical transfer function of the telescope pupil only (psfr.nOtf x psfr.nOtf)
otfStat         :: Optical transfer function of the telescope pupil including the static phase psfr.static_map (psfr.nOtf x psfr.nOtf)
Created by      :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 10/04/2019
                      
Change Record:  ::
------------HEADER END----------------
%}

function obj = computeStaticOpticalTransferFunction(obj,varargin)
inputs = inputParser;
inputs.addRequired('obj',@(x) isa(x,'psfReconstruction'));
inputs.addParameter('statModes',{[]},@iscell);
inputs.addParameter('flagStaticMap',true,@islogical);
inputs.parse(obj,varargin{:});

statModes           = inputs.Results.statModes;
obj.flags.staticMap = inputs.Results.flagStaticMap;

%Verify if adding static modes is needed
if ~isempty(statModes{1})
    obj.statModes = statModes{1};
    obj.statCoefs = statModes{2};
else
    obj.statCoefs = 0;
end

% Define the static aberrations
obj.res.static.full_map = 0;
if obj.flags.staticMap
    obj.res.static.full_map =  obj.trs.tel.pupil.*(obj.trs.tel.static_fitting  ... %+ obj.trs.cam.diff_field_stat
        + obj.res.static.stat_map_interp)*1e-9*2*pi/obj.trs.cam.wavelength;
end

% Adding additional aberrations
if ~isempty(statModes{1})
    nRes = sqrt(size(obj.statModes,1));
    if nRes ~= obj.trs.tel.resolution 
        fprintf('The size of the input modes does not comply with the telescope resolution\n');
    else
        obj.res.static.full_map = obj.res.static.full_map + 1e-9*2*pi/obj.trs.cam.wavelength*obj.trs.tel.pupil.*reshape(obj.statModes*obj.statCoefs',nRes,[]);
    end
end

% Get the OTF
P               = puakoTools.enlargePupil(obj.trs.tel.pupil,2);
phi_stat        = puakoTools.enlargePupil(obj.res.static.full_map ,2);
E               = P.*exp(1i*phi_stat);
otfStat         = fftshift(puakoTools.fftCorrel(E,E));
obj.otf.otfStat = otfStat /max(otfStat(:));
