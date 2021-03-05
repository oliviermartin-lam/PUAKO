
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NGS CASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par files
parFileKeckScaoStructPuako
parm.sci.x = 0;% arcsec
parm.cam.exposureTime = 1000;

% launch the AO simulation
aoSys = aoSystem(parm,'runSimulation',true);

% instantiate the telemetry 
trs = telemetry(aoSys,'none','none',struct(),{''},'flagNoisemethod','nonoise');
r0_estFromTelem = trs.res.seeing.r0;
L0_estFromTelem = trs.res.seeing.L0;


%% INSTANTIATE A PSFR CLASS

% hack the estimation to check to set aside PSFR innacuracy on account of
% incorrect atmospheric parameter estimation
psfr = psfReconstruction(trs,'flagAoPattern','square','flagNoisemethod','nonoise');

% diplay results: comparing the sytem PSF and the reconstructed PSF
displayResults(psfr)

errorBreakDown(psfr)

% psfr = psfReconstruction(trs,'flagAoPattern','square','flagDphiMethod','zonal');
% 
% % diplay results: comparing the sytem PSF and the reconstructed PSF
% displayResults(psfr)
% 
% errorBreakDown(psfr)


%% INSTANTIATE A PSFR CLASS

% hack the estimation to check to set aside PSFR innacuracy on account of
% incorrect atmospheric parameter estimation
trs.res.seeing.r0 = parm.atm.r0;
trs.res.seeing.L0 = parm.atm.L0;

psfrBench = psfReconstruction(trs,'flagAoPattern','square','flagNoiseMethod','nonoise');

% diplay results: comparing the sytem PSF and the reconstructed PSF
displayResults(psfrBench)

errorBreakDown(psfrBench)


%% PRIME FOR THE HYBRID RECONSTRUCTION (ESTIMATION OF 
pr = prime(psfr,'fitCn2',true);
displayResults(pr)
errorBreakDown(pr)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LGS CASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par files
parFileKeckSlaoStructPuako
%parFileSLAO
parm.cam.exposureTime = 1000;
%parm.sci.x = 0;% arcsec

% launch the AO simulation
aoSys = aoSystem(parm,'runSimulation',true);

% instantiate the telemetry 
trs = telemetry(aoSys,'none','none',struct(),{''},'flagNoiseMethod','nonoise');
r0_estFromTelem = trs.res.seeing.r0;
L0_estFromTelem = trs.res.seeing.L0;

%% INSTANTIATE A PSFR CLASS

psfr = psfReconstruction(trs,'flagAoPattern','square','flagNoiseMethod','nonoise');

% diplay results: comparing the sytem PSF and the reconstructed PSF
displayResults(psfr)

errorBreakDown(psfr)

%% PRIME FOR THE HYBRID RECONSTRUCTION (ESTIMATION OF 
pr = prime(psfr);
displayResults(pr)
errorBreakDown(pr)


%% INSTANTIATE A PSFR CLASS

% hack the estimation to check to set aside PSFR innacuracy on account of
% incorrect atmospheric parameter estimation
trs.res.seeing.r0 = parm.atm.r0;
trs.res.seeing.L0 = parm.atm.L0;

psfrBench = psfReconstruction(trs,'flagAoPattern','square');

% diplay results: comparing the sytem PSF and the reconstructed PSF
displayResults(psfrBench)

errorBreakDown(psfrBench)


