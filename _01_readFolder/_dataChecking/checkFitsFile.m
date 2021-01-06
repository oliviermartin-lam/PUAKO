function [out,flagStatus] = checkFitsFile(obj,file_fits)
inputs = inputParser;
inputs.addRequired('obj', @(x) isa(x,'puako'));
inputs.addRequired('file_fits', @(x) ischar(x) | isstruct(x) | iscell(x));
inputs.parse(obj,file_fits);

if isempty(obj.data_folder_path)
    uiwait(errordlg('You must provide an existing folder path first by setting the path_imag property before selecting a .fits file'));
    flagStatus = 'error';
    out = [];
    return;
end

% Instantiating
flagStatus = 'ok';
if ~iscell(file_fits)
    file_fits = {file_fits};
end
nIMAG = numel(file_fits);
flagStatusFITS = cell(1,nIMAG);
out    = cell(1,nIMAG);
nC     = zeros(1,nIMAG);

% Check if the .fits file exist
folder_fits = dir(obj.data_folder_path.imag);

for k=1:nIMAG
    idx = find(contains({folder_fits.name},file_fits{k}));
    isfits = any(idx);
    if isfits
        nC(k) = numel({folder_fits(idx).name});
        tmp = cell(1,nC(k) );
        for j=1:nC(k) 
            tmp{j} =  [obj.path_imag,folder_fits(idx(j)).name];            
        end
        out{k} = tmp;
        if numel(idx) >= 1
            fprintf(['Selection of ',num2str(numel(idx)),' .fits file','s'*(numel(idx) > 1),'\n']);
            flagStatusFITS{k} = 'ok';
        else
            flagStatusFITS{k} = 'warning';
        end
    else
        flagStatusFITS{k} = 'warning';
    end
end


if any(strcmp(flagStatusFITS,'warning'))
    warning(sprintf('%d .fits file have not been selected because they do not exist.',nnz(strcmp(flagStatusFITS,'warning'))));
        flagStatus = 'warning';
end

% Put results into a single cell
nF = sum(nC);
tmp = cell(1,nF);
j=1;
for k=1:nIMAG
    tmp(1,j:j+nC(k)-1) = out{k};
    j = j+nC(k);
end
out = tmp;