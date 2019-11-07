%{
------------HEADER-----------------
Objective          ::  Select the valid actuators regarding the outer and
inner pupils size
INPUT VARS
dmModes         :: the acutators commands
dmCom            :: Outer and inner pupil diameter in meter
tel                :: DM actuators pitch
Cn                    :: Number of 1D actuators over the full pupil x/y-axis
fitL0
best
initR0
D1
D2
badModesList
aoMode
nMin
nMax

OUTPUT VARS
r0            :: The estimated r0
L0            :: The estimated L0
dr0            :: The 3-sigma precision on the r0
dL0            :: The 3-sigma precision on the L0
validZernikeModes   :: The j-index vector of Zernikes modes used for fitting
z_coefs_mod            :: The corresponding fitted coefficients
z_coefs_var            :: The empirical Zernike modes vaariance
(denoising/dealiasing and removed from bad modes)

Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 11/01/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function varargout = getr0L0FromDMcommands(dmModes,dmCom,tel,Cn,varargin)
inputs = inputParser;
inputs.addRequired('dmModes',@isnumeric);
inputs.addRequired('dmCom',@isnumeric );
inputs.addRequired('tel',@isstruct);
inputs.addRequired('Cn',@isnumeric);
inputs.addParameter('fitL0',true,@islogical);
inputs.addParameter('flagBest',false,@islogical);
inputs.addParameter('flagMedian',false,@islogical);
inputs.addParameter('initR0',[],@isnumeric);
inputs.addParameter('D1',11.25,@isnumeric);
inputs.addParameter('D2',0,@isnumeric);
inputs.addParameter('badModesList',[],@isnumeric);
inputs.addParameter('aoMode','NGS',@ischar);
inputs.addParameter('nMin',4,@isnumeric);
inputs.addParameter('nMax',120,@isnumeric);
inputs.parse(dmModes,dmCom,tel,Cn,varargin{:});

%1\ Parsing inputs
fitL0    = inputs.Results.fitL0;
flagBest = inputs.Results.flagBest;
flagMedian = inputs.Results.flagMedian;
initR0  = inputs.Results.initR0;
D1       = inputs.Results.D1;
D2       = inputs.Results.D2;
nMin     = inputs.Results.nMin;
nMax     = inputs.Results.nMax;

%2\ Estimating the Zernike modes rms values
%2.1 Defining the Zernike modes over the real telescope pupil
nRes   = sqrt(size(dmModes,1));
validZernikeModes = nMin:nMax;
zern_tel  = zernike(validZernikeModes,nRes);
z_modes = zern_tel.modes;

% %2.2 Defining the Zernike modes over the truncated pupil
% if D1 ~=0
%     % truncated pupil definition    
%     res = getGridCoordinates(nRes,nRes,tel.Ddm/2);   
%     P         = double(logical(tools.interpolate(tel.pupil,nRes)));
%     P_trunc = P.*(res.r2D<=D1/2).*(res.r2D>D2/2);
%     % define Zernike modes over the truncated pupil
%     zern_trunc   = zernike(validZernikeModes,nRes,'pupil',P_trunc);    
%     % Calculate the projection matrix
%     zTrunc_modes= zern_trunc.modes;
%     proj   = pinv(zTrunc_modes'*zTrunc_modes)*zTrunc_modes'*z_modes;
% else
%     proj = 1;
%     D1 = tel.Dcircle;
% end

%3\ Process the telemetry to estimate the Zernike modes
%3.1 Actuators commands to Zernike
u2z      = pinv(z_modes)*dmModes;
z_coefs = u2z*dmCom;
%3.2. Denoising
z_noise_var = diag(u2z*Cn*u2z');
z_coefs_var   = std(z_coefs,[],2).^2 - z_noise_var;
%3.3. Dealiasing
%polynomial model of the variance excess due to the aliasing
pp = [0.99893514401310313 7.7066823200766521e-4,...
    -1.8233364176012401e-4 1.9037601227012146e-5 -2.7620673576089771e-7];
[nz,mz] = nmOrder(validZernikeModes);
aliasExcess = 0;
for k=1:numel(pp)
    aliasExcess = aliasExcess+pp(k).*nz.^(k-1);
end
z_coefs_var = z_coefs_var./aliasExcess';

%3.4. Bad modes list
if isempty(inputs.Results.badModesList)
    nValid = numel(validZernikeModes);
    badModesList = ones(1,nValid);
    if strcmp(inputs.Results.aoMode,'NGS')
        wbm = [[4:6],[87,88 101,102 74,75, 62,63, 51 41 33 25 116 117 115]];
    else
        wbm = [[4:6],[7,8,11,16,17,18,22,29,30,37,38,...
            44,46,47,56,60,62,64,66,67,68,79,81,82,84,86,...
            88,90,92,93,94,95,106,110,112,114,116]];
    end
    badModesList(intersect(validZernikeModes,wbm)) = 0;
    badModesList(mz <=1) = 0;
    badModesList = logical(badModesList);
end

% Removing bad modes
validZernikeModes = validZernikeModes(badModesList);
z_coefs_var = z_coefs_var(badModesList);
z_noise_var = z_noise_var(badModesList);

%4\ Model-fitting the Zernike modes
if flagBest % Multiple fitting and keep the most precise one    
    % instantiation
    nM   = nMin:nMax/4; iN   = length(nM);
    r0   = zeros(1,iN);dr0  = zeros(1,iN);
    L0   = zeros(1,iN);dL0  = zeros(1,iN);
    % loop over
    for k=1:iN
        jindex = validZernikeModes(k:end);
        [r0(k),L0(k),dr0(k),dL0(k)] = fitZernikeVariance(z_coefs_var(k:end),jindex,D1,'fitL0',fitL0,'initR0',initR0);
    end    
    % Measure the joint precision over adjusted parameters
    eps   = hypot(dr0./r0,dL0./L0);
    iBest = find(eps == min(eps));
    z_noise_var = z_noise_var(iBest:end);
    z_coefs_var(iBest:end) = z_coefs_var(iBest:end);
    % Keep the most precise result
    r0    = r0(iBest);dr0   = dr0(iBest);
    L0    = L0(iBest); dL0   = dL0(iBest);
    
elseif flagMedian % Multiple fitting with median removal
    
    k=0;
    jindex = validZernikeModes;
    zz = z_coefs_var;
    eps = 1;
    while eps > 1e-1
        k = k+1;
        [r0,L0,dr0,dL0] = fitZernikeVariance(zz,jindex,D1,'fitL0',fitL0,'initR0',initR0);
        z_coefs_mod = zernikeVarianceModel([D1/r0,D1/L0],jindex);
        z_diff = abs((zz - z_coefs_mod)./zz);
        eps = std(z_diff);
        idx_good = abs(z_diff) < 3*eps;
        jindex = jindex(idx_good);
        zz = zz(idx_good);
        
        if k >=50 || numel(jindex) < 50
            break;
        end
    end
    z_noise_var = z_noise_var(idx_good);
    z_coefs_var = z_coefs_var(idx_good);
    
else    % Single fitting over the validZernikeModes modes
    [r0,L0,dr0,dL0] = fitZernikeVariance(z_coefs_var,validZernikeModes,D1,'fitL0',fitL0,'initR0',initR0);
    jindex = validZernikeModes;
end

%5\ Gathering outputs
varargout{1} = r0;
varargout{2} = L0;
varargout{3} = dr0;
varargout{4} = dL0;
varargout{5} = jindex';
varargout{6} = zernikeVarianceModel([D1/r0,D1/L0],jindex);
varargout{7} =z_coefs_var;
varargout{8} =z_noise_var;
