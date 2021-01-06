function flagStatus = readTRSFolder(obj,path_trs)
inputs = inputParser;
inputs.addRequired('obj', @(x) isa(x,'puako'));
inputs.addRequired('path_trs', @ischar);
inputs.parse(obj,path_trs);

flagStatus = 'ok';

%% 1\ Check if at least one of the folder exist
if ~isfolder(path_trs)
    uiwait(errordlg('You must provide an existing folder path'));
    flagStatus = 'error';
    return;
end
%% 2\ Keep in memory

%2.2. TRS folder
if isfolder(path_trs)
    obj.data_folder_path.trs = path_trs;
    % Check emptiness of folders
    folder_trs = dir(obj.data_folder_path.trs);
    issav = any(find(contains({folder_trs.name},'.sav')));
    if ~issav
        warning('There''re no .sav file in the TRS directory, please check it out.');
        flagStatus = 'warning';
    end
else
    warning(['The folder',path_trs,' does not exist\n']);
    flagStatus = 'warning';
end





