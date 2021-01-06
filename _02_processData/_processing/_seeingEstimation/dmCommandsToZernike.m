function [u2z] = dmCommandsToZernike(dmModes,jMin,jMax,varargin)
inputs = inputParser;
inputs.addRequired('dmModes',@isnumeric);
inputs.addRequired('Jmin',@isnumeric);
inputs.addRequired('Jmax',@isnumeric);
inputs.addParameter('mskPup',true(size(dmModes,1)),@islogical);
inputs.addParameter('mskModes',true(size(dmModes,2)),@islogical);
inputs.parse(dmModes,jMin,jMax,varargin{:});

mskPup = inputs.Results.mskPup;
mskModes = inputs.Results.mskModes;

% Zernike modes defined on the circular pupil
nRes = sqrt(size(dmModes,1));
zDM   = zernike_puako(jMin:jMax,nRes);
zDM = zDM.modes;

nRes_trunc = sqrt(length(mskPup(mskPup)));
zTrunc   = zernike_puako(jMin:jMax,nRes_trunc);
zD1  = 0*zDM;
zD1(mskPup(:),:) = zTrunc.modes;

%u2z = pinv(zD1'*zD1)*zD1'*zDM*pinv(zDM)*dmModes; %
u2z = pinv(zD1'*zD1)*zD1'*dmModes; %
u2z = u2z(:,mskModes(:));
return

%Pupil mask
% if inputs.Results.Dpup == 0
%     Dpup = 10.54*pitch;
% end    
% res = getGridCoordinates(nRes,nRes,Dpup/2);
% msk = (res.r2D<=D1/2).*(res.r2D>D2/2);
% 
% %DM modes spatial covariance matrix
% matII = dmModes'*dmModes/nRes;
% 
% % Zernike modes defined on the circular pupil and then masked
% z   = zernike_puako(jMin:jMax,nRes,'pupil',msk);
% z = z.modes;
% 
% % New zernike modes
% appZernike = dmModes*pinv(full(matII))'*(z'*dmModes)'/nRes;
% 
% %Projection
% matZZ = appZernike'*appZernike/nRes;
% matZI = appZernike'*dmModes/nRes;
% u2z = pinv(matZZ,inputs.Results.tol)*matZI;
