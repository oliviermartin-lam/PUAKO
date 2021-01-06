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
elseif contains(msg,'omartin')
    path_workspace = '/home/omartin/Projects/KVS/CODES/PUAKO/';
    addpath(genpath(path_workspace));
    path_root = '/run/media/omartin/8C8C-B02C/KECK_DATA/';
    path_calib = [path_root,'CALIBRATION/'];
    path_save='/home/omartin/Projects/KVS/_results/seeingComparison/';
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
p.getPrimePSF('path_save',path_save,'fitBg',true,'fitGains',false);

%% Display results;
res = statisticalAnalysis(path_save);
res.displayPlot('TRSR0','PRIR0','legend1','Telemetry $r_0$ (m)','legend2','Image $r_0$ (m)','xyline',true)