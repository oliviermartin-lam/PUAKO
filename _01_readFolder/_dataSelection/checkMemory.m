function flagStatus  = checkMemory(obj,id_trs,id_imag)
inputs = inputParser;
inputs.addRequired('obj',@(x) isa(x,'puako'));
inputs.addRequired('id_trs',@isstruct);
inputs.addRequired('id_imag',@isstruct);
inputs.parse(obj,id_trs,id_imag);

switch obj.userChoice
    case {1|2|5}
        dataSize = sum([id_trs.size]) + sum([id_imag.size]);
    case {3|4}
        dataSize = sum([id_imag.size]);
end
if isunix
    [~,msg] = unix('free -m --si');
    msg = split(msg);
    sysMem = str2num(msg{14});
    rat = 1e2*dataSize/sysMem;
    flagStatus = 'ok';
    if rat < 40
        fprintf('I''ll use %.1f %s of the available system memory\n',rat,'%');
    elseif rat>=40 && rat<70
        fprintf('Caution: I''m going to use %.1 %s of the available system memory\n',rat,'%');
    elseif rat>=70 && rat < 95
        x =  input(['Warning: I need ',num2str(rat),'% of the available system memory to load all the data, do you want me to continue (''y''/''n'') ?']);
        flagStatus = 'error';
        if isstring(x) || ischar(x)
            if contains(upper(x),'YES') || contains(upper(x),'Y')
                flagStatus = 'ok';
            end
        elseif isnumeric(x) && x
            flagStatus = 'ok';
        end
    else
        uiwait(errordlg(['You have only ',num2str(sysMem/1e3,2),'Gb of available system memory and I need ' num2str(dataSize/1e3,2),...
            'Gb to load all the data, I abort the loading process. Please check data_id.trs.size to select a sub-sample of files.']));
        return
    end
else
    warning(sprintf(['\n You''re not using a UNIX plateforme so the data file size estimation is compromised.\n']));
    flagStatus = 'warning';
end
