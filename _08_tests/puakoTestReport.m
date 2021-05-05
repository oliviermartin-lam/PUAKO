
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NGS CASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par files
parFileKeckScaoStructPuako
parm.sci.x = 0;% arcsec
parm.cam.exposureTime = 5000;

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

%% INSTANTIATE A PSFR CLASS

% hack the estimation to check to set aside PSFR innacuracy on account of
% incorrect atmospheric parameter estimation
airmass = 1/cos(parm.atm.zenithAngle);
trs.res.seeing.r0 = parm.atm.r0 * airmass^(-3/5);
trs.res.seeing.L0 = parm.atm.L0;
psfrBench = psfReconstruction(trs,'flagAoPattern','circle','flagNoiseMethod','nonoise');

% diplay results: comparing the sytem PSF and the reconstructed PSF
displayResults(psfrBench)

errorBreakDown(psfrBench,'display',true);


%% PRIME FOR THE HYBRID RECONSTRUCTION (ESTIMATION OF 
pr = prime(psfrBench,'fitCn2',false,'fitBg',false);
displayResults(pr)
errorBreakDown(pr,'display',true);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LGS CASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par files
parFileKeckSlaoStructPuako
parm.cam.exposureTime = 5000;
parm.sci.x = 0;% arcsec
parm.nGs.x = 60;% arcsec

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

%% PRIME FOR THE HYBRID RECONSTRUCTION (ESTIMATION OF 
pr = prime(psfrBench,'fitCn2',false,'fitBg',false);
displayResults(pr)
errorBreakDown(pr,'display',true);


