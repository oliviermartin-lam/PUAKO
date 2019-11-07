function [n_tot,idxgood,flagStatus] = getDataNumberFromUserSpecification(dataID,id)
inputs = inputParser;
inputs.addRequired('dataID', @isstruct);
inputs.addRequired('id', @(x) iscell(x) |isnumeric(x) | ischar(x) | isstruct(x));
inputs.parse(dataID,id);

%1\ Check whether the folder has been read already
flagStatus = 'ok';
n_tot = 0;
if isempty(dataID)
    uiwait(errordlg('You must load data first using the readFolder function.'));
    flagStatus = 'error';
    return
end
if ~iscell(id)
    id = {id};
end

%2\ Get the number of files the user want to load
  n = numel(id);      
  n_tot = 0;
  pos_tmp = [];
  obj_name_tmp = 'none';
  idxgood = zeros(1,n);
    for j=1:n
        id_tmp = id{j};       
        if isnumeric(id_tmp)
            obj_name = dataID.imag(id_tmp).name;
            if ~all(contains(obj_name,obj_name_tmp))
                n_tot = n_tot + sum(~contains(obj_name,obj_name_tmp));
                obj_name_tmp = obj_name;
                idxgood(j) = 1;
            end
        elseif ischar(id_tmp) || isstruct(id_tmp)            
            if isempty(pos_tmp)  || (~isempty(pos_tmp)  && any(pos_tmp ~= (contains({dataID.trs.id},id_tmp) & contains({dataID.imag.id},id_tmp))) )% avoid redundant call
                pos_tmp = contains({dataID.trs.id},id_tmp) & contains({dataID.imag.id},id_tmp);
                n_tot = n_tot + sum(pos_tmp);
                obj_name_tmp = dataID.imag(pos_tmp).name;
                idxgood(j) = 1;
            end
        elseif iscell(id_tmp)
            [n_tmp,flagStatus] = getDataNumberFromUserSpecification(dataID,id_tmp);
            n_tot = n_tot + n_tmp;
        else
            uiwait(errordlg('The format is not recognized.'));
            data_index = [];         
            flagStatus = 'error';
            return            
        end                
    end
  
    