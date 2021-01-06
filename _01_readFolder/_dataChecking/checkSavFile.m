function [out,flagStatus] = checkSavFile(obj,file_sav)
inputs = inputParser;
inputs.addRequired('obj', @(x) isa(x,'puako'));
inputs.addRequired('file_sav', @(x) ischar(x) | isstruct(x) | iscell(x));
inputs.parse(obj,file_sav);

if isempty(obj.data_folder_path)
    uiwait(errordlg('You must provide an existing folder path first by setting the path_trs property before selecting a .sav file'));
    flagStatus = 'error';
    out = [];
    return;
end

flagStatus = 'ok';

if ~iscell(file_sav)
    file_sav = {file_sav};
end
nTRS = numel(file_sav);
out    = cell(1,nTRS);
flagStatusSAV = cell(1,nTRS);
nC     = zeros(1,nTRS);


% Check if the .sav file exist
folder_sav = dir(obj.data_folder_path.trs);
for k=1:nTRS
    idx = find(contains({folder_sav.name},file_sav{k}));
    issav = any(idx);
    
    if issav
        nC(k) = numel({folder_sav(idx).name});
        tmp = cell(1,nC(k) );
        for j=1:nC(k) 
            tmp{j} =  [obj.path_trs,folder_sav(idx(j)).name];            
        end
        out{k} = tmp;        
        if numel(idx) >= 1
            fprintf(['Selection of ',num2str(numel(idx)),' .sav file','s'*(numel(idx) > 1),'\n']);
             flagStatusSAV{k} = 'ok';
        else
             flagStatusSAV{k} = 'warning';
        end
    end
end
if any(strcmp(flagStatusSAV,'warning'))
    warning(sprintf('%d .sav file have not been selected because they do not exist.',nnz(strcmp(flagStatusSAV,'warning'))));
        flagStatus = 'warning';
end

% Put results into a single cell
nF = sum(nC);
tmp = cell(1,nF);
j=1;
for k=1:nTRS
    tmp(1,j:j+nC(k)-1) = out{k};
    j = j+nC(k);
end
out = tmp;