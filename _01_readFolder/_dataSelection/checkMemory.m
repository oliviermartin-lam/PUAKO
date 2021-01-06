function flagStatus  = checkMemory(dataSize)
inputs = inputParser;
inputs.addRequired('dataSize',@isnumeric);
inputs.parse(dataSize);

if isunix
    [~,msg] = unix('free -m --si');
    msg = split(msg);
    sysMem = str2num(msg{14});
    rat = ceil(1e2*dataSize/sysMem);
    flagStatus = 'ok';
    if rat < 40
        fprintf('I''ll use %.1f %s of the available system memory\n',rat,'%');
    elseif rat>=40 && rat<70
        fprintf('Caution: I''m going to use %.1f %s of the available system memory\n',rat,'%');
    elseif rat>=70 && rat < 95
        x =  input(['Warning: I need ',num2str(rat,1),'% of the available system memory to load all the data, do you want me to continue (''y''/''n'') ?']);
        if contains(upper(x),'Y')
            flagStatus = 'ok';
        else
            fprintf('Ok, so I abort the loading\n');
            flagStatus = 'error';
            return
        end
    else
        uiwait(errordlg(['You have only ',num2str(sysMem/1e3,2),'Gb of available system memory and I need ' num2str(dataSize/1e3,2),...
            'Gb to load all the data, I abort the loading process. Please check file_fits and file_sav properties to select a sub-sample of files.']));
        flagStatus = 'error';
        return
    end
else
    warning(sprintf(['\n You''re not using a UNIX plateforme so the data file size estimation is compromised.\n']));
    flagStatus = 'warning';
end
