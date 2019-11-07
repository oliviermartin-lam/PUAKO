function [id_trs,id_imag,fitsHdr,flagStatus] = selectFilesFromUser(obj,id)
if ~exist('id')
    id = 'ALL';
end
inputs = inputParser;
inputs.addRequired('obj',@(x) isa(x,'puako'));
inputs.addRequired('id',@(x) isstring(id) | ischar(id) | isnumeric(id) | iscell(id) );
inputs.parse(obj,id);

% Determine which data should be restored: .sav, .fits or both
flagStatus = 'ok';
if isempty(obj.data_id)
    uiwait(errordlg('You must load data first using the readFolder function.'));
    flagStatus = 'error';
    return
else
    dataID_trs = obj.data_id.trs;
    dataID_imag = obj.data_id.imag;
    if isempty([dataID_trs.id])
        dataID_trs.id = dataID_imag.id;
    end
    if isempty([dataID_imag.id])
        dataID_imag.id = dataID_trs;
    end
end

%% id is numeric
% Grab the corresponding data ID
if ((ischar(id) || isstring(id)) && strcmpi(id,'ALL')) ||(iscell(id) && contains(upper(id),'ALL'))
    % Process all data in obj.data_id.trs.id
    id_trs = dataID_trs;
    id_imag = dataID_imag;
    fitsHdr = obj.fitsHeader;
elseif (ischar(id) || isstring(id)) && ~strcmpi(id,'ALL')
    % Process a specific file by its name (id = 'n0004')
    obj_pos = contains({dataID_trs.id},id) & contains({dataID_imag.id},id) ;
    if ~isempty(obj_pos)
        id_trs = dataID_trs(obj_pos);
        id_imag = dataID_imag(obj_pos);
        fitsHdr = obj.fitsHeader(obj_pos);
        if contains({dataID_trs.id},id) ~= contains({dataID_imag.id},id)
            warning(sprintf(['\nThere are not as many IMAG than TRS files.\n']));
            flagStatus = 'warning';
        end
    else
        uiwait(errordlg('The file name you''ve specified does not exist or has not been loaded.'));
        flagStatus = 'error';
        return
    end
    
    %% id is numeric
elseif isnumeric(id)
    % Process a specific file by its id (id = 5)
    id_trs = dataID_trs(id);
    id_imag = dataID_imag(id);
    fitsHdr = obj.fitsHeader(id);
    
    %% id is a cell
elseif iscell(id)
    % Determine all the data files covered by the user's request
    [n_tot,idxgood,flagStatus] = getDataNumberFromUserSpecification(obj.data_id,id);
    if strcmp(flagStatus,'error')
        flagStatus = 'error';
        return
    end    
    % Fill in the id structures
    id_trs = struct('id',char(n_tot),'aomode',char(n_tot),'size',zeros(1,n_tot),'name',char(n_tot));
    id_imag = struct('id',char(n_tot),'aomode',char(n_tot),'size',zeros(1,n_tot),'name',char(n_tot));
    fitsHdr =struct('field',cell(1,n_tot),'value',cell(1,n_tot),'comments',cell(1,n_tot));        
    
    %Looping on files
    for k=1:nnz(idxgood)
        id_tmp = id{idxgood(k)};
        if isnumeric(id_tmp)              
            id_trs(k+nnz(obj_pos)-1) = dataID_trs(id_tmp);
            id_imag(k+nnz(obj_pos)-1) = dataID_imag(id_tmp);
            fitsHdr(k+nnz(obj_pos)-1) = obj.fitsHeader(id_tmp);
        elseif (ischar(id_tmp) || isstring(id_tmp))
            obj_pos = contains({dataID_trs.id},id_tmp) & contains({dataID_imag.id},id_tmp) ;
            if ~isempty(obj_pos)
                id_trs(k:k+nnz(obj_pos)-1) = dataID_trs(obj_pos);
                id_imag(k:k+nnz(obj_pos)-1) = dataID_imag(obj_pos);
                fitsHdr(k:k+nnz(obj_pos)-1) = obj.fitsHeader(obj_pos);
                if any(contains({dataID_trs.id},id_tmp) ~= contains({dataID_imag.id},id_tmp))
                    warning(sprintf(['\nThere are not as many IMAG than TRS files.\n']));
                    flagStatus = 'warning';
                end
            end
        else
            uiwait(errordlg('The file name you''ve specified does not exist or has not been loaded.'));
            flagStatus = 'error';
            return
        end
    end
    %% if format is not recognized
else
    uiwait(errordlg('You must provide a valid file identity like: ''ALL'', 1, [4:10] or ''n0004''.'));
    flagStatus = 'error';
    return
end

nObj = numel({id_trs.id});
if nObj == 1
    fprintf('I''ll load 1 data file then.\n');
else
    fprintf('I''ll load %d data files then.\n',nObj);
end
