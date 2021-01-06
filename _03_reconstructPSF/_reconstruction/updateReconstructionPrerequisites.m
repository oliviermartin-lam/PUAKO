function psfr = updateReconstructionPrerequisites(psfr,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfReconstruction'));
inputs.parse(psfr,varargin{:});

%1\ Update the sampling
[psfr.otf.dk,psfr.otf.nTimes,psfr.trs.tel.resolution,psfr.otf.nOtf] = defineSampling(psfr.trs.dm.nActuators,psfr.trs.cam.samp,psfr.otf.fov_fit);
            
%2\ Update the instrumental features size
if ~isempty(psfr.trs.tel.static_map(:))
    psfr.res.static.stat_map_interp = puakoTools.interpolate(psfr.trs.tel.static_map,psfr.trs.tel.resolution);
    psfr.trs.tel.pupil              = double(logical(psfr.res.static.stat_map_interp));
else
    psfr.res.static.stat_map_interp = zeros(psfr.trs.tel.resolution);
    psfr.trs.tel.pupil              = puakoTools.interpolate(psfr.trs.tel.pupil,psfr.trs.tel.resolution,'nearest');
end
if ~isempty(psfr.trs.cam.diff_field_stat(:))
    psfr.trs.cam.diff_field_stat    = puakoTools.interpolate(psfr.trs.cam.diff_field_stat,psfr.trs.tel.resolution);
else
    psfr.trs.cam.diff_field_stat    = zeros(psfr.trs.tel.resolution);
end
if ~isempty(psfr.trs.tel.static_fitting(:))
    psfr.trs.tel.static_fitting     = puakoTools.interpolate(psfr.trs.tel.static_fitting,psfr.trs.tel.resolution);
else
    psfr.trs.tel.static_fitting     = zeros(psfr.trs.tel.resolution);
end

for iSrc = 1:numel(psfr.trs.sky)
    [psfr.trs.sky(iSrc).SR,psfr.trs.sky(iSrc).dSR]  = puakoTools.getStrehl(psfr.trs.cam.image(:,:,iSrc),psfr.trs.tel.pupil,psfr.trs.cam.samp);
end

% Get the Pupil OTF
P               = puakoTools.enlargePupil(psfr.trs.tel.pupil ,2);
psfr.otf.otfDL  = fftshift(puakoTools.fftCorrel(P,P));
psfr.otf.otfDL  = psfr.otf.otfDL /max(psfr.otf.otfDL(:));

%3\ Update the influence function
psfr = updateInfluenceFunction(psfr);