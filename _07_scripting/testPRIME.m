clear all;close all;clc;
date = '20130801';

[~,msg] = unix('echo "$USER"');
if contains(msg,'psfr')    
    path_workspace = '/home/psfr/PUAKO/PUAKO/';
    addpath(genpath(path_workspace));
    path_root = '/net/vm-psfr/psfrData/NIRC2/';
    path_calib = '/net/vm-psfr/psfrData/CALIBRATION/';
    path_save='/net/vm-psfr/psfrData/PUAKO_results/';
    path_imag = [path_root,date,'/IMAG/'];   
    path_trs = [path_root,date,'/TRS/'];   
    path_ncpa = [];
elseif contains(msg,'omartin')
    path_workspace = '/home/omartin/Projects/KVS/CODES/PUAKO/';
    addpath(genpath(path_workspace));
    path_root = '/run/media/omartin/OlivierMartinHDD/DATA/KECK_DATA/';
    path_calib = [path_root,'CALIBRATION/'];
    path_save='/home/omartin/Projects/KVS/_results/psfComparison/';
    path_imag = [path_root,date,'/IMAG/'];   
    path_trs = [path_root,date,'/TRS/'];   
    path_save = [path_save,date,'/'];
    path_ncpa = [path_calib,'AOSYSTEM/STATIC/phase_diversity_results_20200311_average_PD3.sav'];
elseif contains(msg,'sragland')
    path_workspace = '/home/sragland/PUAKO/PUAKO/';
    addpath(genpath(path_workspace));
    path_root = '/net/vm-psfr/psfrData/NIRC2/';
    path_calib = '/net/vm-psfr/psfrData/CALIBRATION/';
    path_save='/net/vm-psfr/psfrData/PUAKO_RESULTS_SAM/RESULTS/';    
    path_imag = [path_root,date,'/IMAG/'];   
    path_trs = [path_root,date,'/TRS/'];   
    path_ncpa = [];
end

p = puako('path_imag',path_imag,'path_trs',path_trs,'path_calibration',path_calib);

%% Get PSFR
close all;

%p.getPrimePSF('objname',{'n0059'},'resolution',150,'fov',110,'path_ncpa',path_ncpa);
p.getPrimePSF('objname',{'n0004'},'resolution',150,'fov',170,'ron',300);%,'fitStatModes',[1:36],'statModesFunction','piston');
p.vizualizeResults('resultsof','prime')
1e2*p.psfp.catalog_fit.fvu

%p.getPrimePSF('path_save',path_save,'fitBg',true,'fitGains',[true true true]);



