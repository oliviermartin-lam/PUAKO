date = '20170314';

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
    path_root = '/run/media/omartin/8C8C-B02C/KECK_DATA/';
    path_calib = [path_root,'CALIBRATION/'];
    path_save='/home/omartin/Projects/KVS/_results/psfComparison/';
    path_imag = [path_root,date,'/IMAG/'];   
    path_trs = [path_root,date,'/TRS/'];   
    path_save = [path_save,date,'/'];
elseif contains(msg,'sragland')
    path_workspace = '/home/sragland/PUAKO/PUAKO/';
    addpath(genpath(path_workspace));
    path_root = '/net/vm-psfr/psfrData/NIRC2/';
    path_calib = '/net/vm-psfr/psfrData/CALIBRATION/';
    path_save='/net/vm-psfr/psfrData/PUAKO_RESULTS_SAM/RESULTS/';    
    path_imag = [path_root,date,'/IMAG/'];   
    path_trs = [path_root,date,'/TRS/'];   
end
%% Get PSFR
p = puako('path_imag',path_imag,'path_trs',path_trs,'path_calibration',path_calib);
p.getRecPSF('path_save',path_save,'fitBg',true);

