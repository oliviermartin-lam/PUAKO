function selectSamples(obj,objname,varargin)
% Check if the folder has been read already
if isempty(obj.data_folder_path)
    uiwait(errordlg('You must provide an existing folder path first by setting the path_trs property before processing .sav files'));
    obj.flagStatus = 'error';
    return;
end
% Select all files if no inputs
if ~exist('objname','var')
    objname = {'all'};
end

inputs = inputParser;
inputs.addRequired('obj',@(x) isa(x,'puako'));
inputs.addRequired('objname', @(x) iscell(x) | ischar(x));
inputs.addParameter('path_imag', obj.data_folder_path.imag, @ischar);
inputs.addParameter('path_trs', obj.data_folder_path.trs, @ischar);
inputs.addParameter('path_calibration',obj.data_folder_path.calibration, @ischar);
inputs.addParameter('path_save',[], @ischar);
inputs.addParameter('yesToAll', true, @islogical);
inputs.addParameter('threshSteeringMirror', 20, @isnumeric);
inputs.addParameter('getImageOnly', false, @islogical);

inputs.parse(obj,objname,varargin{:});
objname = inputs.Results.objname;
getImageOnly = inputs.Results.getImageOnly;

%1\ Checking data folder if different from the initial config
if ~strcmp(obj.path_imag,obj.data_folder_path.imag)
    obj.path_imag =  inputs.Results.path_imag;
end
if ~strcmp(obj.path_trs,obj.data_folder_path.trs)
    obj.path_trs =  inputs.Results.path_trs;
end
if ~strcmp(obj.path_calibration,obj.data_folder_path.calibration)
    obj.path_calibration =  inputs.Results.path_calibration;
end

%2\ Sub-select samples
if ~strcmpi(objname,'ALL')
    obj.file_fits = {objname};
    obj.file_sav = {objname};
else
    obj.file_fits = '.fits';
    obj.file_sav = '.sav';
end

% 3\ Check the size in Mb
if obj.keepInMemory
    dataSize_imag = getSampleSize(obj.data_folder_id.imag,obj.file_fits);
    obj.flagStatus = checkMemory(dataSize_imag);
    if ~getImageOnly
        dataSize_trs = getSampleSize(obj.data_folder_id.trs,obj.file_sav);
        obj.flagStatus = checkMemory(dataSize_imag+ dataSize_trs);
    end
    if strcmp(obj.flagStatus,'error')
        return
    end
end

%% TAKE THE OVERLAPPING AREA IF USER WANTS TRS + IMAG DATA
nTRS = numel(obj.file_sav);
    nIMAG = numel(obj.file_fits);
    
    trsname = cell(1,nTRS);
    for k=1:nTRS
        tmp = split(obj.file_sav{k},'/');
        tmp = split(tmp{end},'.');
        tmp = split(tmp{1},'_');
        trsname{k} = tmp{1};
    end
    
    imagname = cell(1,nIMAG);
    for k=1:nIMAG
        tmp = split(obj.file_fits{k},'/');
        tmp = split(tmp{end},'.');
        tmp = split(tmp{1},'_');
        imagname{k} = tmp{1};
    end
    
if ~getImageOnly        
    if nTRS~=0 && nTRS~=nIMAG
        fprintf('Sub-selection of the TRS/IMAG pairs\n');
        name_inter = intersect({obj.data_folder_id.trs.name},imagname);
        sub_trs = find(contains({obj.data_folder_id.trs.name},name_inter));
        sub_imag= find(contains({obj.data_folder_id.imag.name},name_inter));
        trsname = {obj.data_folder_id.trs(sub_trs).name};
        imagname = {obj.data_folder_id.imag(sub_imag).name};
        selectSamples(obj,trsname);
    end
    
    obj.trs_idx = find(contains({obj.data_folder_id.trs.name},trsname));
    obj.imag_idx = find(contains({obj.data_folder_id.imag.name},imagname));
else
    obj.trs_idx = [];
    obj.imag_idx = find(contains({obj.data_folder_id.imag.name},imagname));
end
            
