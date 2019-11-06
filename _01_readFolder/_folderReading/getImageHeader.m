function flagStatus = getImageHeader(obj,varargin)
inputs = inputParser;
inputs.addRequired('obj', @(x) isa(x,'puako'));
inputs.parse(obj,varargin{:});

% Read the fits header
nIMAG = numel(obj.data_id.imag);
hdrfield = cell(1,nIMAG);
hdrvalue = cell(1,nIMAG);
hdrcomments = cell(1,nIMAG);
warnFlag = 0;

tic();
if nIMAG == 1
    fprintf('Reading the .fits file header...')
    tmp = [obj.folder_data.imag];
    
    if isfile(tmp) % check if the file exist
        tmp = fitsinfo(tmp);
        if isfield(tmp,'PrimaryData')
            hdr = tmp.PrimaryData.Keywords;
            hdrfield{1} = hdr(:,1);
            hdrvalue{1} = hdr(:,2);
            hdrcomments{1} = hdr(:,3);
        else
            fprintf('\n');
            warning(sprintf('\nLooks like the file %s format is ackward and there is no PrimaryData field, sorry !',tmp));
            flagStatus = 'warning';
            warnFlag = 1;
        end
    else
        fprintf('\n');
        warning(sprintf('\nThe file %s does not exist, please check',tmp));
        flagStatus = 'warning';
        warnFlag = 1;
    end
    
else
    fprintf('Reading %d .fits files headers...',nIMAG);
    
     
            
    for k=1:nIMAG
        tmp = [obj.folder_data.imag,obj.data_id.imag(k).id];
        if isfile(tmp) % check if the file exist
            tmp = fitsinfo(tmp);
            if isfield(tmp,'PrimaryData')
                hdr = tmp.PrimaryData.Keywords;
                hdrfield{k} = hdr(:,1);
                hdrvalue{k} = hdr(:,2);
                hdrcomments{k} = hdr(:,3);
            else
                warnFlag = warnFlag + 1;
                fprintf('\n');
                warning(sprintf('\nLooks like the file %s format is ackward and there is no PrimaryData field, sorry !',tmp));
                flagStatus = 'warning';
            end
        else
            fprintf('\n');
            warning(sprintf('\nThe file %s does not exist, please check',tmp));
            flagStatus = 'warning';
            warnFlag = warnFlag + 1;
        end
    end
end

t0 = toc();
if warnFlag == 1
    fprintf('...Done in %.2g s. You''ve got 1 warning, please check what''s wrong.\n',t0)
elseif warnFlag > 1
    fprintf('...Done in %.2g s. You''ve got %d warnings, please check what''s wrong.\n',t0,warnFlag)
else
    fprintf('...Done in %.2g s. You''ve got no warnings, you rock buddy.\n',t0)
    flagStatus = 'ok';
end

obj.fitsHeader = struct('field',hdrfield,'value',hdrvalue,'comments',hdrcomments);
end
