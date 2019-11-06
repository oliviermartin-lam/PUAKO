function flagStatus = getDataID(obj,varargin)
inputs = inputParser;
inputs.addRequired('obj', @(x) isa(x,'puako'));
inputs.parse(obj,varargin{:});

% Get the files name
nTRS = 0;
if isfield(obj.folder_data,'trs') && ~isempty(obj.folder_data.trs)
    trs_list = dir(obj.folder_data.trs);
    nTRS = numel(trs_list);
end
nIMAG = 0;
if isfield(obj.folder_data,'imag') && ~isempty(obj.folder_data.imag)
    imag_list = dir(obj.folder_data.imag);
    nIMAG = numel(imag_list);
end

if nTRS
    if nIMAG>nTRS
        warning(sprintf(['\nThere are more IMAG data than TRS files !\n']));
        flagStatus = 'warning';
    elseif nIMAG<nTRS
        warning(sprintf(['\nThere are more TRS data than IMAG files !\n']));
        flagStatus = 'warning';
    end
end

% Fill in the data_id structures
if nTRS > 1
    %remove the ./ and../ folders
    trs_list = trs_list(3:end);
    nTRS = numel(trs_list);
    %Get the AO observing mode
    aomode = cell(1,nTRS);
    trsBytes = cell(1,nTRS);
    imagBytes = cell(1,nTRS);
    obj_name  = cell(1,nTRS);
    
    for k=1:nTRS
        trsBytes{k} = double(trs_list(k).bytes)/1000^2;
        imagBytes{k} = double(imag_list(k).bytes)/1000^2;
        aomode{k} = 'NGS';
        if isempty(strfind(trs_list(k).name,'NGS'))
            aomode{k} = 'LGS';                                
        end        
        tmp = split(trs_list(k).name,'_');
        obj_name{k} = tmp{1};
    end
    
    % Store TRS data ID
    obj.data_id.trs  = struct('id',{trs_list.name},'aomode',aomode,'size',trsBytes,'name',obj_name);
    flagStatus = 'ok';

elseif nTRS == 1
    aomode = 'NGS';
    if ~isempty(strfind(obj.folder_data.trs,'LGS'))
        aomode = 'LGS';
    end
    obj.data_id.trs.id = obj.folder_data.trs;
    obj.data_id.trs.aomode = aomode;
    obj.data_id.trs.size =  double(trs_list.bytes)/1000^2;
    tmp = split(trs_list.name,'_');
    tmp = split(tmp{1},'.');
    obj.data_id.trs.name =  tmp{1};
    flagStatus = 'ok';
else
    obj.data_id.trs.id = [];
    obj.data_id.trs.aomode = [];
    obj.data_id.trs.size = [];
    obj.data_id.trs.name = [];
    flagStatus = 'ok';
end


if nIMAG > 1
    imag_list = imag_list(3:end);
    nIMAG = numel(imag_list);
    imagBytes = cell(1,nIMAG);
    for k=1:nIMAG
        imagBytes{k} = double(imag_list(k).bytes)/1000^2;
        tmp = split(imag_list(k).name,'_');
        tmp = split(tmp{1},'.');
        obj_name{k} = tmp{1};
    end
    aomode = cell(1,nIMAG);
    if nTRS 
        obj.data_id.imag  = struct('id',{imag_list.name},'aomode',aomode,'size',imagBytes,'name',obj_name);
    else
        % do not know the aomode without the TRS file
        obj.data_id.imag  = struct('id',{imag_list.name},'size',imagBytes,'name',obj_name);
    end
    flagStatus = 'ok';

elseif nIMAG == 1
    obj.data_id.imag.id = obj.folder_data.imag;
    if nTRS
        obj.data_id.imag.aomode = obj.data_id.trs.aomode;
    else
        obj.data_id.imag.aomode = 'unknown';
    end
    obj.data_id.imag.size =  double(imag_list.bytes)/1000^2;
    tmp = split(imag_list.name,'_');
    tmp = split(tmp{1},'.');
    obj.data_id.imag.name = tmp{1};
    flagStatus = 'ok';

else
    obj.data_id.imag = [];
    warning(sprintf(['\nThere is IMAG data loaded !\n']));
    flagStatus = 'warning';
end
