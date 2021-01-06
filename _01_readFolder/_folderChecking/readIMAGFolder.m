function flagStatus = readIMAGFolder(obj,path_imag,varargin)
inputs = inputParser;
inputs.addRequired('obj', @(x) isa(x,'puako'));
inputs.addRequired('path_imag', @ischar);
inputs.parse(obj,path_imag,varargin{:});

flagStatus = 'ok';
%% 1\ Check if at least one of the folder exist
if ~isfolder(path_imag)
    uiwait(errordlg('You must provide an existing folder path'));
    flagStatus = 'error';
    return;    
end

%% 2\ Keep in memory

if isfolder(path_imag)
    obj.data_folder_path.imag = path_imag;
    % Check emptiness of folders
      folder_imag = dir(obj.data_folder_path.imag);
    isfits = any(find(contains({folder_imag.name},'.fits')));
    if ~isfits
        warning(['There''re no .fits file in the IMAG directory, please check it out.']);
        flagStatus = 'warning';
    end
else
    obj.data_folder_id.imag = [];
    warning(['The folder',path_imag,' does not exist']);
end

%% Create DARK and SKY folders
if ~obj.noToAll
    flagStatus = sortDarkandSkyandDivFiles(path_imag,'threshSteeringMirror',obj.threshSteeringMirror,'yesToAll',obj.yesToAll);
end