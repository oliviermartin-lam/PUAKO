function flagStatus = getTRSDataID(obj,varargin)
inputs = inputParser;
inputs.addRequired('obj', @(x) isa(x,'puako'));
inputs.addParameter('path_trs', obj.path_trs,@ischar);
inputs.parse(obj,varargin{:});
flagStatus = 'ok';

% Check the folder
path_tel = inputs.Results.path_trs;

% Read files inside
trs_list = dir(path_tel);

% Get .fits file only
trs_list = trs_list(contains({trs_list.name},'.sav'));
if isempty(trs_list)
    uiwait(errordlg('You must provide an existing folder path'));
    flagStatus = 'error';
end
nTRS = numel(trs_list);

% Count the number of .fits file
path_trs = cell(1,nTRS);
name_trs  = cell(1,nTRS);
type_trs = cell(1,nTRS);
imagBytes = cell(1,nTRS);
hdr = cell(1,nTRS);
aomode = cell(1,nTRS);
flagStatusHDR= cell(1,nTRS);
for k=1:nTRS
    % Get the path
    path_trs{k} = [path_tel,trs_list(k).name];
    % Get the name, type and ao mode
    tmp = split(trs_list(k).name,'_');
    name_trs{k} = tmp{1};
    % Get the AO mode
    aomode{k} = 'NGS';
    if isempty(strfind(tmp{2},'NGS'))
        aomode{k} = 'LGS';
    end
    tmp = split(tmp{3},'.');
    type_trs{k} = tmp{2};
    % Get the size in Mega Bytes
    imagBytes{k} = double(trs_list(k).bytes)/1000^2;
    
    % Get the header from data_folder_imag
    flagStatusHDR{k} = 'warning';
    if ~isempty(obj.data_folder_id.imag)
        kGood = strcmp({obj.data_folder_id.imag.origin_name},name_trs{k});
        if any(kGood)
            hdr{k} = obj.data_folder_id.imag(find(kGood)).hdr;
            flagStatusHDR{k} = 'ok';
        else
            if ~isempty(obj.data_folder_id.dark)
                % Maybe the file has been moved to a Dark  directory
                kGood = strcmp({obj.data_folder_id.dark.origin_name},name_trs{k});
                if any(kGood)
                    hdr{k} = obj.data_folder_id.dark(find(kGood)).hdr;
                     flagStatusHDR{k} = 'ok';
                else
                    if ~isempty(obj.data_folder_id.sky)
                        % Maybe the file has been moved to a Dark  directory
                        kGood = strcmp({obj.data_folder_id.sky.origin_name},name_trs{k});
                        if any(kGood)
                            hdr{k} = obj.data_folder_id.sky(find(kGood)).hdr;                        
                             flagStatusHDR{k} = 'ok';
                        end
                    else
                        if ~isempty(obj.data_folder_id.div)
                            % Maybe the file has been moved to a Dark  directory
                            kGood = strcmp({obj.data_folder_id.div.origin_name},name_trs{k});
                            if any(kGood)
                                hdr{k} = obj.data_folder_id.div(find(kGood)).hdr;
                                flagStatusHDR{k} = 'ok';
                            end
                        end
                    end
                end
            end
        end                               
    end
end

if any(strcmp( flagStatusHDR,'warning'))
    warning(sprintf('%d TRS files are not associated with a fits header !',nnz(strcmp(flagStatusHDR,'warning'))));
    flagStatus = 'warning';
end
        
obj.data_folder_id.trs  = struct('name',name_trs,'type',type_trs,'aomode',aomode,'size',imagBytes,'hdr',hdr,'path',path_trs);
