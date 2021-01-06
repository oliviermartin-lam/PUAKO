function [dataSize,flagStatus]  = getSampleSize(path_id,path_file)
inputs = inputParser;
inputs.addRequired('path_id',@isstruct);
inputs.addRequired('path_file',@iscell);
inputs.parse(path_id,path_file);

%1\ Get the size of selected files
idxgood = contains({path_id.path},path_file);
if any(idxgood)
    dataSize = sum([path_id(idxgood).size]);
else
    dataSize = 0;
    warning('No files selected');
    flagStatus = 'warning';
    return;
end