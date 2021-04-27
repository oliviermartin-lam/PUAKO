
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NGS CASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par files
parFileKeckScaoStructPuako
parm.sci.x = 10;% arcsec
parm.cam.exposureTime = 100;

% launch the AO simulation
aoSys = aoSystem(parm,'runSimulation',true);

% instantiate the telemetry 
trs = telemetry(aoSys,'none','none',struct(),{''},'flagNoisemethod','nonoise');
r0_estFromTelem = trs.res.seeing.r0;
L0_estFromTelem = trs.res.seeing.L0;


%% INSTANTIATE A PSFR CLASS

% hack the estimation to check to set aside PSFR innacuracy on account of
% incorrect atmospheric parameter estimation
psfr = psfReconstruction(trs,'flagAoPattern','square','flagNoisemethod','nonoise','flagResidualMethod','slopes-based');

% diplay results: comparing the sytem PSF and the reconstructed PSF
displayResults(psfr)

errorBreakDown(psfr,'display',true);

% psfr = psfReconstruction(trs,'flagAoPattern','square','flagDphiMethod','zonal');
% 
% % diplay results: comparing the sytem PSF and the reconstructed PSF
% displayResults(psfr)
% 
% errorBreakDown(psfr)


% case on-axis
%Wavelength		2.18 micron
%Strehl Image		72% +/- 0.086	
%Strehl Marechal 	74.1%	
%Strehl Parenti		74.2%	
%Wavefront error (nm)	190	
%-------------------------------
%Residual Static		0
%Atmospheric Fitting	141.2
%-------------------------------
%Servo-lag		66.4
%WFS noise		0
%WFS aliasing		80.13
%-------------------------------
%Tip-tilt bandwidth	73.19
%Tip-tilt noise		0
%-------------------------------
%Total Anisoplanatism	0
%Focal anisoplanatism	0
%Angular-anisoplanatism	0
%Anisokinetism		0
%-------------------------------

% CASE 10"-off-axis

%Wavelength		2.18 micron
%Strehl Image		13.9% +/- 0.028	
%Strehl Marechal 	12.3%	
%Strehl Parenti		12.7%	
%Wavefront error (nm)	502.3	
%-------------------------------
%Residual Static		0
%Atmospheric Fitting	138.4
%-------------------------------
%Servo-lag		71
%WFS noise		0
%WFS aliasing		78.55
%-------------------------------
%Tip-tilt bandwidth	71.58
%Tip-tilt noise		0
%-------------------------------
%Total Anisoplanatism	465.6
%Focal anisoplanatism	0
%Angular-anisoplanatism	465.6
%Anisokinetism		0


%% INSTANTIATE A PSFR CLASS

% hack the estimation to check to set aside PSFR innacuracy on account of
% incorrect atmospheric parameter estimation
airmass = 1/cos(parm.atm.zenithAngle);
trs.res.seeing.r0 = parm.atm.r0 * airmass^(-3/5);
trs.res.seeing.L0 = parm.atm.L0;

psfrBench = psfReconstruction(trs,'flagAoPattern','square','flagNoiseMethod','nonoise');

% diplay results: comparing the sytem PSF and the reconstructed PSF
displayResults(psfrBench)

errorBreakDown(psfrBench)


%% PRIME FOR THE HYBRID RECONSTRUCTION (ESTIMATION OF 
pr = prime(psfr,'fitCn2',false);
displayResults(pr)
errorBreakDown(pr,'display',true);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LGS CASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par files
parFileKeckSlaoStructPuako
parm.cam.exposureTime = 1000;
parm.sci.x = 10;% arcsec
parm.nGs.x = 0;% arcsec

% launch the AO simulation
aoSys = aoSystem(parm,'runSimulation',true);

% instantiate the telemetry 
trs = telemetry(aoSys,'none','none',struct(),{''},'flagNoiseMethod','nonoise');
r0_estFromTelem = trs.res.seeing.r0;
L0_estFromTelem = trs.res.seeing.L0;

%% INSTANTIATE A PSFR CLASS

psfr = psfReconstruction(trs,'flagAoPattern','square','flagNoiseMethod','nonoise','flagAnisoMethod','flicker');

% diplay results: comparing the sytem PSF and the reconstructed PSF
displayResults(psfr)

errorBreakDown(psfr,'display',true);

%% PRIME FOR THE HYBRID RECONSTRUCTION (ESTIMATION OF 
pr = prime(psfr);
displayResults(pr)
errorBreakDown(pr)


%% INSTANTIATE A PSFR CLASS

% hack the estimation to check to set aside PSFR innacuracy on account of
% incorrect atmospheric parameter estimation
airmass = 1/cos(parm.atm.zenithAngle);
trs.res.seeing.r0 = parm.atm.r0 * airmass^(-3/5);
trs.res.seeing.L0 = parm.atm.L0;
psfrBench = psfReconstruction(trs,'flagAoPattern','square','flagNoiseMethod','nonoise');

% diplay results: comparing the sytem PSF and the reconstructed PSF
displayResults(psfrBench)

errorBreakDown(psfrBench,'display',true);


