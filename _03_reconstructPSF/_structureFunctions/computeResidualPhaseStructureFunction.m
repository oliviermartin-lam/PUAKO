%{
------------HEADER-----------------
Objective           ::  Compute bi-dimensional phase structure function of the tip-tilt excluded residual phase

INPUT VARS
psfr                :: The psfReconstruction top-level class
jZernGain           :: (optionnal) list of Zernike modes (ex: 4:6 for focus + astigmatisms terms) to enable modal gains retrieval
OUTPUT VARS
sf_2D               :: bi-dimensional phase structure function map (Toeplitz) of the aliasing error (psfr.nOtf x psfr.nOtf)
Cho                 :: Covariance matrix of the tip-tilt excluded residual error in the acturators space (nActu^2 x nActu^2)
Created by          :: O. Beltramo-Martin - ONERA/LAM
Creation date       :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function obj = computeResidualPhaseStructureFunction(obj,varargin)
inputs = inputParser;
inputs.addRequired('obj',@(x) isa(x,'psfReconstruction'));
inputs.addParameter('jZernGain',[],@isnumeric);
inputs.addParameter('method','Vii',@ischar);
inputs.addParameter('flagResidualMethod','slopes-based',@ischar);
inputs.parse(obj,varargin{:});
obj.trs.mat.jZernGain = inputs.Results.jZernGain;

%1\ Get the covariance matrix in the actuators space
if strcmp(inputs.Results.flagResidualMethod,'dm-based')
    du = diff(obj.trs.dm.com,[],2)/obj.trs.holoop.gain;
    obj.cov.Cao  = du*du'/obj.trs.wfs.nExp;
else
    obj.cov.Cao  = obj.trs.rec.res*obj.trs.rec.res'/obj.trs.wfs.nExp;    
end

%2\ Get the phase structure function in the actuators space
Du     = (2*pi/obj.trs.cam.wavelength)^2*(obj.cov.Cao  - obj.trs.res.noise(1).Cn_ho);

%3\ Get the point-wise phase structure function assuming stationarity
if strcmp(inputs.Results.method,'zonal_hr')
    [obj.otf.otfAO,obj.sf.Dao]= modes2Otf(Du,obj.trs.mat.dmIF_hr,obj.trs.tel.pupil,obj.otf.nOtf,'method','zonal');
elseif strcmp(inputs.Results.method,'zonal')
        [obj.otf.otfAO,obj.sf.Dao]= modes2Otf(Du,obj.trs.mat.dmIF,obj.trs.tel.pupil,obj.otf.nOtf,'method','zonal');
else
    [~,obj.sf.Dao]= modes2Otf(Du,obj.trs.mat.dmIF_hr,obj.trs.tel.pupil,obj.otf.nOtf,'method','Vii');
end


%4\ Get modal phase structure functions
%if the user provides a list of Zernike modes, this code will calculate the residual phase structure function for each mode +
%    the remaining one

if ~isempty(obj.trs.mat.jZernGain)
    %4.1 Defining the modal basis, here Zernike
    nZ = numel(obj.trs.mat.jZernGain);
    z  = zernike_puako(obj.trs.mat.jZernGain,obj.trs.dm.nActuators);
    %4.2 Instantiating
    sfFit_2D_full           = obj.sf.Dao;
    n                       = size(obj.sf.Dao,1);
    obj.sf.Dao_z            = zeros(n,n,nZ+1);
    obj.sf.Dao_z(:,:,end)   = sfFit_2D_full;
    obj.cov.Cao_z           = zeros(obj.trs.dm.nActuators^2,obj.trs.dm.nActuators^2,nZ+1);
    obj.trs.mat.Hz          = cell(1,nZ);
    %4.3 Get the DM spatial filter
    bif     = xineticsInfluenceFunctionModel(obj.trs.dm.pitch); % specific to Keck
    bif.setInfluenceFunction(obj.trs.dm.nActuators,obj.trs.dm.nActuators);
    dmFilter= pinv(full(bif.modes),1e-1);
    
    %4.4 loop on considered modes.
    for k=1:nZ
        %4.4.1 Calculate the filtering matrix
        zz                      = dmFilter*squeeze(z.modes(:,k));
        obj.trs.mat.Hz{k}       = zz*pinv(zz); % keep only the kth mode
        %4.4.2 Get the phase structure functions in the actuators space
        obj.cov.Cao_z(:,:,k)    = obj.trs.mat.Hz{k}*(obj.cov.Cao - obj.trs.res.noise.Cn_ho)*obj.trs.mat.Hz{k}';
        Cu                      = (2*pi/obj.trs.cam.wavelength)^2*obj.cov.Cao_z(:,:,k) ;
        %4.4.3 Get the point-wise phase structure funciton assuming stationarity
        [~,obj.sf.Dao_z(:,:,k)] = modes2Otf(Cu,obj.trs.mat.dmIF_hr,obj.trs.tel.pupil,obj.otf.nOtf,'method','Vii');
    end
    %4.5 Get the remaining residual SF subtracted from zernike modes.
    obj.sf.Dao_z(:,:,k+1)   = obj.sf.Dao - squeeze(sum(obj.sf.Dao_z(:,:,1:end-1),3));
    obj.sf.Cao_z(:,:,k+1)   = obj.cov.Cao - squeeze(sum(obj.cov.Cao_z(:,:,1:end-1),3));
else
    obj.sf.Dao_z    = obj.sf.Dao;
    obj.cov.Cao_z   = obj.cov.Cao;
end

