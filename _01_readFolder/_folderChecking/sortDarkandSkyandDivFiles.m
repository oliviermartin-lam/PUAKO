function flagStatus = sortDarkandSkyandDivFiles(path_image,varargin)
inputs = inputParser;
inputs.addRequired('path_image', @ischar);
inputs.addParameter('threshSteeringMirror', 20,@isnumeric);
inputs.addParameter('yesToAll', true,@islogical);
inputs.parse(path_image,varargin{:});
threshSteeringMirror = inputs.Results.threshSteeringMirror;
yesToAll = inputs.Results.yesToAll;

if ~isfolder(path_image) && ~isfile(path_image)
    fprintf('Something went wrong with the data id.\n');
    flagStatus = 'error';
    return;
end

if isfile(path_image) && ~isfolder(path_image)
    fprintf('I do not sort the data folder.\n');
    flagStatus = 'warning';
    return;
end

if ~isunix
    fprintf('Sorry, your plateform is not supported');
    flagStatus = 'warning';
    return
else
    % Create a Dark folder
    if ~isfolder([path_image,'DARK'])  && ~isfolder([path_image,'Dark']) && ~isfolder([path_image,'dark'])
        msg = input('There is no Dark folder within the .fits folder, do you want me to create it ?  (''y''/''n'')');
        if contains(upper(msg),'Y')
            unix(['mkdir ',path_image,'Dark/']);            
        end
    end
    
    % Create a Sky folder
    if ~isfolder([path_image,'SKY'])  && ~isfolder([path_image,'Sky']) && ~isfolder([path_image,'sky'])
        msg = input('There is no Sky folder within the .fits folder, do you want me to create it ?  (''y''/''n'')');
        if contains(upper(msg),'Y')
            unix(['mkdir ',path_image,'Sky/']);
        end
    end
    
    % Create a Div folder
    if ~isfolder([path_image,'DIV'])  && ~isfolder([path_image,'Div']) && ~isfolder([path_image,'div'])
        msg = input('There is no Diversity folder within the .fits folder, do you want me to create it ?  (''y''/''n'')');
        if contains(upper(msg),'Y')
            unix(['mkdir ',path_image,'Div/']);
        end
    end
end


%% Move dark files in the dar directory
msg = 'Y';
if ~yesToAll
    msg = input('Do you want me to sort dark/sky/diversity files ?  (''y''/''n'')');
end
if contains(upper(msg),'Y')
    if ~yesToAll
        msg = input('Ok, and can I change the filename wih the format type_itime_coadds_date_n ? Files with same itime/coadds an date won''t be erased but incremented  (''y''/''n'')');
    end
    if contains(upper(msg),'Y')
        changeNamePermission = true;
    else
        changeNamePermission = false;
    end
    
    img_name = dir(path_image);
    fits_name = {img_name(contains({img_name.name},'.fits')).name};
    nFits = numel(fits_name);
    dark_name_tmp = cell(1,nFits);
    sky_name_tmp = cell(1,nFits);
    div_name_tmp = cell(1,nFits);
    jBg = 1;
    jSky = 1;
    jDiv = 1;
    nDark = 0;
    nSky = 0;
    nDiv = 0;
    
    % loop on fits file
    for k=1:nFits
        fits_path = [path_image,fits_name{k}];
        hdr = fitsinfo(fits_path);
        hdr = hdr.PrimaryData.Keywords;
        % Read the keywords
        shutter = upper(cell2mat(hdr(strcmp(hdr(:,1),'SHRNAME'),2)));
        itime = num2str(cell2mat(hdr(strcmp(hdr(:,1),'ITIME'),2)));
        itime = strrep(itime,'.','p');
        coadds = num2str(cell2mat(hdr(strcmp(hdr(:,1),'COADDS'),2)));
        dateFile = num2str(cell2mat(hdr(strcmp(hdr(:,1),'DATE-OBS'),2)));
        AOFMX = str2num(cell2mat(hdr(strcmp(hdr(:,1),'AOFMX'),2)));
        AOFMY = str2num(cell2mat(hdr(strcmp(hdr(:,1),'AOFMY'),2)));
        OBWF = str2num(cell2mat(hdr(strcmp(hdr(:,1),'OBWF'),2)));
        LSPROP = cell2mat(hdr(strcmp(hdr(:,1),'LSPROP'),2));
        if strcmp(LSPROP,'yes') 
            AOFCLGFO = str2num(cell2mat(hdr(strcmp(hdr(:,1),'AOFCLGFO'),2)))*1e3; %mm)
            wfs_dz = abs(OBWF - AOFCLGFO); % not sure it works
        else
            AOFCNGFO = str2num(cell2mat(hdr(strcmp(hdr(:,1),'AOFCNGFO'),2)))*1e3; %mm)
            wfs_dz = abs(OBWF - AOFCNGFO);
        end
        
        
        % 1\ Check for DARK: the detector shutter is closed
        if contains(shutter,'CLOSE')
            nDark = nDark + 1;
            % this is a dark
            dark_name = ['dark_itime_',itime,'_coadds_',coadds,'_date_',dateFile,'.fits'];
            if any(strcmp(dark_name,dark_name_tmp))
                good = 0;
                while good
                    jBg = jBg+1;
                    dark_name = ['dark_itime_',itime,'_coadds_',coadds,'_date_',dateFile,'_',num2str(jBg),'.fits'];
                    if ~isfile([path_image,'Dark/',dark_name]) % check if the files does not alreadye exist
                        good = 1;
                    end
                end
            end
            dark_name_tmp{k} = dark_name;
            if changeNamePermission
                unix(['mv ',path_image,fits_name{k},' ',path_image,'Dark/',dark_name]);
            else
                unix(['mv ',path_image,fits_name{k},' ',path_image,'Dark/',fits_name{k}]);
            end
        end
        
        % 2\ Check for SKY BACKGROUND: the sterring mirror is largely displaced        
        if hypot(AOFMX,AOFMY) > threshSteeringMirror %supposely the sterring miror was used to get a sky
            nSky = nSky + 1;
            % this is a sky
            sky_name = ['sky_itime_',itime,'_coadds_',coadds,'_date_',dateFile,'.fits'];
            if any(strcmp(sky_name,sky_name_tmp))
                good = 0;
                while good
                    jSky = jSky+1;
                    sky_name = ['sky_itime_',itime,'_coadds_',coadds,'_date_',dateFile,'_',num2str(jSky),'.fits'];
                    if ~isfile([path_image,'Sky/',dark_name]) % check if the files does not already exist
                        good = 1;
                    end
                end
            end
            sky_name_tmp{k} = sky_name;
            if changeNamePermission
                unix(['mv ',path_image,fits_name{k},' ',path_image,'Sky/',sky_name]);
            else
                unix(['mv ',path_image,fits_name{k},' ',path_image,'Sky/',fits_name{k}]);
            end
        end
        
        
        % 3\ Check for Phase diversity measurements: the WFS stage is move in z-position : AOWFC0
        if wfs_dz > 3  %difference between the WFS stage z-position and the model offset value = f(EL,AZ)
            nDiv  =nDiv + 1;
            % this is a sky
            div_name = ['div_zpos',num2str(wfs_dz),'_itime_',itime,'_coadds_',coadds,'_date_',dateFile,'.fits'];
            if any(strcmp(div_name,div_name_tmp))
                good = 0;
                while good
                    jDiv = jDiv+1;
                    div_name = ['div_zpos',num2str(OBWF - AOFCNGFO),'_itime_',itime,'_coadds_',coadds,'_date_',dateFile,'_',num2str(jDiv),'.fits'];
                    if ~isfile([path_image,'Div/',div_name]) % check if the files does not alreadye exist
                        good = 1;
                    end
                end
            end
            div_name_tmp{k} = div_name;
            if changeNamePermission
                unix(['mv ',path_image,fits_name{k},' ',path_image,'Div/',div_name]);
            else
                unix(['mv ',path_image,fits_name{k},' ',path_image,'Div/',fits_name{k}]);
            end
        end
    end
    if changeNamePermission
        fprintf(['I have moved ',num2str(nDark),' dark, ',num2str(nSky),' sky and ',num2str(nDiv), ' diversity calibration files into corresponding folder with name changing\n']);
    else
        fprintf(['I have moved ',num2str(nDark),' dark, ',num2str(nSky),' sky and ',num2str(nDiv), ' diversity calibration files into corresponding folder without name changing\n']);
    end
end
flagStatus = 'ok';