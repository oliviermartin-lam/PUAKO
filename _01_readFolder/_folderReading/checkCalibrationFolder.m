function obj = checkCalibrationFolder(obj,path_calibration,varargin)
inputs = inputParser;
inputs.addRequired('obj', @(x) isa(x,'puako'));
inputs.addRequired('path_calibration', @ischar);
inputs.parse(obj,path_calibration,varargin{:});

% Grab the CALIBRATION folder path
obj.folder_data.calibration = path_calibration;
if ~ isfolder(obj.folder_data.calibration)
    warning('The calibration path you''ve provided does not exist, you won''t be able to use all my capabilities ! ')
end
% Create the MASSDIMM sub-folder
obj.folder_data.massdimm = [path_calibration,'_massdimm/'];
if ~ isfolder(obj.folder_data.massdimm)
    fprintf('Creating the _massdimm/ sub-folder in the calibration folder.\n')
    mkdir(obj.folder_data.massdimm);
end

