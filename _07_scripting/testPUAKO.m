clear all;
close all;

%% DEFINE DATA PATH
path_root = '/psfrData/NIRC2/';
path_calib = '/psfrData/PUAKO_RESULTS/CALIBRATION/';
path_save='/psfrData/PUAKO_RESULTS/RESULTS/';
path_night = [path_root,'20130801/'];
testPSFR = true;
p = puako();

%% Test data 
if testData
    p = p.grabData('id',{'n0004'},'path_night',path_night,'path_calibration',path_calib,'path_save',[path_save,'20130801/']);
end

%% Test PSFR
if testPSFR
    t0 = tic();
    p.getRecPSF('id',{'all'},'path_night',path_night,'path_calibration',path_calib,'path_save',[path_save,'20130801/']);
    t0 = toc(t0);    
    p.vizualizeResults('resultsof','psfr');
end