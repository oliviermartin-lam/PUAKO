clear all;close all;
path_workspace = '/home/omartin/Projects/KVS/CODES/PUAKO/';
addpath(genpath(path_workspace));
path_root = '/run/media/omartin/8C8C-B02C/KECK_DATA/';
path_calib = [path_root,'CALIBRATION/'];
path_save_root='/home/omartin/Projects/KVS/_results/psfComparison/';
date = {'20130203' '20130801' '20130914' '20170314'};

for k=1:numel(date)
    path_imag = [path_root,date{k},'/IMAG/'];
    path_trs = [path_root,date{k},'/TRS/'];
    path_save = [path_save_root,date{k},'/'];
    
    % Get PSFR
    p = puako('path_imag',path_imag,'path_trs',path_trs,'path_calibration',path_calib,'yesToAll',true);
    p.getPrimePSF('path_save',path_save,'fitBg',true,'fitGains',true,'fov',300);
end

