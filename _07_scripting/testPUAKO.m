clear all;
close all;

%% DEFINE DATA PATH
[~,msg] = unix('echo "$USER"');
if contains(msg,'psfr')    
    path_workspace = '/home/psfr/PUAKO/PUAKO/';
    addpath(genpath(path_workspace));
    path_root = '/net/vm-psfr/psfrData/NIRC2/';
    path_calib = '/net/vm-psfr/psfrData/CALIBRATION/';
    path_save='/net/vm-psfr/psfrData/PUAKO_results/';
    path_imag = [path_root,'20130801/IMAG/'];   
    path_trs = [path_root,'20130801/TRS/'];   
elseif contains(msg,'omartin')
    path_workspace = '/home/omartin/Projects/KVS/CODES/PUAKO/';
    addpath(genpath(path_workspace));
    path_root = '/run/media/omartin/8C8C-B02C/KECK_DATA/';
    path_calib = [path_root,'CALIBRATION/'];
    path_save='/home/omartin/Projects/KVS/_results/';
    path_imag = [path_root,'20130801/IMAG/'];   
    path_trs = [path_root,'20130801/TRS/'];   
elseif contains(msg,'sragland')
    path_workspace = '/home/sragland/PUAKO/PUAKO/';
    addpath(genpath(path_workspace));
    path_root = '/net/vm-psfr/psfrData/NIRC2/';
    path_calib = '/net/vm-psfr/psfrData/CALIBRATION/';
    path_save='/net/vm-psfr/psfrData/PUAKO_RESULTS_SAM/RESULTS/';    
    path_imag = [path_root,'20130801/IMAG/'];   
    path_trs = [path_root,'20130801/TRS/'];   
end

p = puako('path_imag',path_imag,'path_trs',path_trs,'path_calibration',path_calib);

%% Test data 
%p.grabData('objname',{'n0004'},'path_save',[path_save,'20130801/']);

%% Test PSFR
%if testPSFR
    t0 = tic();
    p.getRecPSF('objname',{'n0004'},'path_save',[path_save,'20130801/']);
    t0 = toc(t0);    
    p.vizualizeResults('resultsof','psfr');
%end

%% Test PRIME
%if testPRIME
    t0 = tic();
    p.getPrimePSF('objname',{'n0004'},'path_save',[path_save,'20130801/']);
    t0 = toc(t0);    
    p.vizualizeResults('resultsof','prime');
%end