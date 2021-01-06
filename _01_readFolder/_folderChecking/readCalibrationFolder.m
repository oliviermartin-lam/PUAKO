function flagStatus = readCalibrationFolder(obj,path_calib)
inputs = inputParser;
inputs.addRequired('obj', @(x) isa(x,'puako'));
inputs.addRequired('path_calib', @ischar);
inputs.parse(obj,path_calib);

flagStatus = 'ok';
%% 1\ Check if at least one of the folder exist
if ~isfolder(path_calib)
    uiwait(errordlg('You must provide an existing folder path for the calibration folder'));
    flagStatus = 'error';
    return;    
end

%% 2\ Keep in memory

if isfolder(path_calib)
    obj.data_folder_path.calibration = path_calib;
    
    %MASSDIMM path
    path_massdimm = [path_calib,'/MASSDIMM/'];
    if  isfolder(path_massdimm)
        obj.data_folder_path.massdimm = path_massdimm;
    else
        msg = 'Y';
        if ~obj.yesToAll
            msg = input(['The MASSDIMM folder does not exist, do you want me to create it ?  (''y''/''n'')']);
        end
        if contains(upper(msg),'Y')
            if isunix
                unix(['mkdir ',path_massdimm]);
            else
                fprintf('Sorry, only UNIX plateforme are supported for folders/files manipulation\n');
                flagStatus = 'warning';
            end
        end
    end
    
    % Check emptiness of folders
    folder_calib = dir(obj.data_folder_path.calibration);
    istel = any(find(contains(upper({folder_calib.name}),'TEL')));
    isao = any(find(contains(upper({folder_calib.name}),'AO')));
    if ~istel
        warning(['There''re no TELESCOPE folder inside the calibration folder, this may compromise PSF-R.\n']);
        flagStatus = 'warning';
    end
    if ~isao
        warning(['There''re no AOSYSTEM folder inside the calibration folder, this may compromise PSF-R.\n']);
        flagStatus = 'warning';
    end    
else
    warning(['The folder',path_calib,' does not exist, you won''t be able to use my capabilities.\n']);
    flagStatus = 'error';
    return;
end





    
end



