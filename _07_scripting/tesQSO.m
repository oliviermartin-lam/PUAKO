clear all;close all;clc;
date = '20170315';

[~,msg] = unix('echo "$USER"');
if contains(msg,'psfr')    
    path_workspace = '/home/psfr/PUAKO/PUAKO/';
    addpath(genpath(path_workspace));
    path_root = '/net/vm-psfr/psfrData/NIRC2/';
    path_calib = '/net/vm-psfr/psfrData/CALIBRATION/';
    path_save='/net/vm-psfr/psfrData/PUAKO_results/';
    path_imag = [path_root,date,'/IMAG/'];   
    path_trs = [path_root,date,'/TRS/'];   
elseif contains(msg,'omartin')
    path_workspace = '/home/omartin/Projects/KVS/CODES/PUAKO/';
    addpath(genpath(path_workspace));
    path_root = '/run/media/omartin/OlivierMartinHDD/DATA/KECK_DATA/';
    path_calib= [path_root,'CALIBRATION/'];
    path_save ='/home/omartin/Projects/KVS/_results/QSO/';
    path_imag = [path_root,date,'/IMAG/QSO/'];   
    path_trs  = [path_root,date,'/TRS/'];   
    path_dark = [path_root,date,'/IMAG/Dark/'];
    path_sky  = [path_root,date,'/IMAG/Sky/'];
elseif contains(msg,'sragland')
    path_workspace = '/home/sragland/PUAKO/PUAKO/';
    addpath(genpath(path_workspace));
    path_root = '/net/vm-psfr/psfrData/NIRC2/';
    path_calib = '/net/vm-psfr/psfrData/CALIBRATION/';
    path_save='/net/vm-psfr/psfrData/PUAKO_RESULTS_SAM/RESULTS/';    
    path_imag = [path_root,date,'/IMAG/'];   
    path_trs = [path_root,date,'/TRS/'];   
end
p = puako('path_imag',path_imag,'path_trs',path_trs,'path_calibration',path_calib,'path_dark',path_dark,'path_sky',path_sky);

objname = {p.data_folder_id.imag.name};
nObj = numel(objname);

%% PRIME on the fits binary

prime_param = zeros(43,nObj);
umax        = 5;

for kObj = 10:10
   
    t0 = tic();

    %1\ Running PRIME on the brightest isolated PSF
    % debug the extraction
    p.getPrimePSF('objname',objname{kObj},'fitGains',[false,false,false],...
        'fitR0',false,'fitBg',false,'resolution',150,'fov',170,'ron',300,'flagBrightestStar',true,...
        'fitStatModes',1:36,'statModesFunction','piston','umax',umax);
    %prime_param(:,kObj) = p.psfp.x_final;
    
end