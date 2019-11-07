function psfr = updateInfluenceFunction(psfr,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfReconstruction'));
inputs.parse(psfr,varargin{:});

% Update the influence function
bif          = xineticsInfluenceFunctionModel(psfr.trs.dm.pitch);
bif.setInfluenceFunction(psfr.trs.dm.nActuators,psfr.trs.tel.resolution);
psfr.trs.mat.dmIF_hr = tools.idlToMatlabIndex(bif.modes,psfr.trs.dm.nActuators,true(psfr.trs.dm.nActuators),'modes');
psfr.trs.mat.dmIF_inv_hr = pinv(full(psfr.trs.mat.dmIF_hr));