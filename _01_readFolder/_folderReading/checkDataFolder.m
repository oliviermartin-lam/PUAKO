function flagStatus = checkDataFolder(obj,path_night,varargin)
inputs = inputParser;
inputs.addRequired('obj', @(x) isa(x,'puako'));
inputs.addRequired('path_night', @ischar);
inputs.parse(obj,path_night,varargin{:});

%1\  Identify user's deepest desires
%The user may want to :
%1: load all TRS/IMAG data from a directory
%2: load all either TRS or IMAG data from a directory only
%3: load the TRS and IMAG file of a single data set
%4: load the IMAG file of a single data set

if ~isfolder(path_night) && ~isfile(path_night)
    uiwait(errordlg('The specified path does not refer to any existing folder or file.'));
    flagStatus = 'error';
    return;    
else
    path_night= replace(path_night,'//','/');
end
             
if contains(path_night,'/TRS')  && ~contains(path_night,'.sav')                       
    obj.userChoice = 2;
    fprintf(['So you want me to read the TRS folder ',path_night,', gottcha.\n']);
elseif contains(path_night,'/IMAG') && ~contains(path_night,'.fits')                     
    obj.userChoice = 3;
    fprintf(['So you want me to read the IMAG folder ',path_night,', gottcha.\n']);
elseif contains(path_night,'.fits')  
    obj.userChoice = 4;
    fprintf(['So you want me to load the IMAG file ',path_night,', gottcha.\n']);
elseif contains(path_night,'.sav')                        
    obj.userChoice = 5;
    fprintf(['So you want me to load the TRS file ',path_night,', gottcha.\n']);
else
    obj.userChoice = 1;
    fprintf(['So you want me to read the full folder ',path_night,', gottcha.\n']);
end

% check data folder accordingly user's choice
switch obj.userChoice
    case {1 ,2}        
        %1\ Grab TRS folder path
        if obj.userChoice == 1
            path_root = path_night;
            obj.folder_data.trs = [path_night,'/TRS/'];
            obj.folder_data.imag = [path_night,'/IMAG/'];
            
            if ~isfolder(obj.folder_data.trs)
                uiwait(errordlg('The TRS folder path does not exist.'));
                flagStatus = 'error';
                return;
            end                  
        else
            path_night = [path_night,'/'];
            path_night = replace(path_night,'//','/');
            obj.folder_data.trs = [path_night];
            tmp = split(path_night,'/');
            path_root = strjoin({tmp{1:end-2}},'/');
            obj.folder_data.imag = [path_root,'/IMAG/'];
        end
        
        % Grab IMAG folder path
        if ~isfolder(obj.folder_data.imag)
            uiwait(errordlg('The IMAG folder path does not exist.'));
            flagStatus = 'error';
            return;
        end
                             
        % Get the date
        path_root= [path_root,'/'];
        path_root = replace(path_root,'//','/');
        tmp = split(path_root,'/');
        obj.folder_data.date =  tmp{end-1};
        obj.folder_data.path_root = path_root;
    
    case 3
        % Grab IMAG folder path
        path_root = path_night;
        obj.folder_data.imag = [path_night,'/'];
        if ~ isfolder(obj.folder_data.imag)
            uiwait(errordlg('The IMAG folder path you''ve specified does not exist.'));
            flagStatus = 'error';
            return;
        end
        % Get the date
        obj.folder_data.path_root = [path_root,'/'];
        tmp = split(path_night,'/');
        obj.folder_data.date =  tmp{end-1};
        
    case 4
        if isfile(path_night)
            % The user has provided the .fits file, no need to get the telemetry
            obj.folder_data.imag = path_night;
            obj.folder_data.trs = [];
        else
            uiwait(errordlg('The file you''ve specified the name does not exist.'));
            flagStatus = 'error';
            return;
        end
        
    case 5
        
        if isfile(path_night)
            % The user has provided the .sav file
            obj.folder_data.trs = path_night;
            % Retrieve the corresponding fits file            
            tmp = split(path_night,'/');
            path_root =  [strjoin({tmp{1:end-2}},'/'),'/'];
            % object name
            tmp = split(tmp{end},'.');
            obj_name = tmp{1};
            if contains(obj_name,'_')
                tmp = split(obj_name,'_');
                obj_name = tmp{1};
            end
            obj.folder_data.imag = [path_root,'IMAG/',obj_name,'.fits'];
            % Get the date
            obj.folder_data.path_root = path_root;
            tmp = split(path_night,'/');
            obj.folder_data.date =  tmp{end-2};
        else
            uiwait(errordlg('The file you''ve specified the name does not exist.'));
            flagStatus = 'error';
            return;
        end
end
flagStatus = 'ok';