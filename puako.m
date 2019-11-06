classdef puako < handle
    
    properties
        % Data ID
        folder_data = [];                                         % Data path structure
        data_id = [];                                               % Data ID structure
        fitsHeader = [];                                           % Fits file header
        % Linked sub-classes
        trs;                                                             % telemetry reading and processing
        psfr;                                                            % Forward psf Reconstruction
        bestFit;                                                       % Hybrid psf reconstruction
        % Results
        res = [];
    end
    
    properties (Hidden = true)
        userChoice;
        readingStatus = false;
    end
    
    methods
        
        function obj = puako()
            
            fprintf('\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
            fprintf('\t--\t\t Hello there, I''m PUAKO and I shall be your guide in this cruel AO world !\t\t\t--\n\n')
            fprintf(['\t--\t You can request me to do whatever is in my capabilities (lucky you), which are: \t\t\t--\n'...
                '\t--\t 1. Data reading\t: load IMAG/TRS files within a night folder and read fits header\t\t--\n'...
                '\t--\t 2. Data processing\t: process telemetry and image data\t\t\t\t--\n'...
                '\t--\t 3. PSF reconstruction\t: reconstruct the PSF in a forward mode (no best-fitting)\t\t\t--\n'...
                '\t--\t 4. Hybrid PSF-R\t: hybrid PSF reconstruction using pupil anf focal plane measurements\t\t--\n'...
                '\t--\t 5. Data analysis\t: display tools for analyzing AO performance and reconstruction results\t\t--'])
            fprintf('\n\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
            
            % Init
            obj.folder_data.trs = [];
            obj.folder_data.imag = [];
            obj.folder_data.date = [];
            obj.folder_data.path_root = [];
            obj.folder_data.calibration = [];
            obj.folder_data.massdimm = [];            
        end
        
        
        %% Data reading
        function readFolder(obj,path_night,path_calibration,varargin)
            inputs = inputParser;
            inputs.addRequired('obj', @(x) isa(x,'puako'));
            inputs.addRequired('path_night', @ischar);
            inputs.addRequired('path_calibration', @ischar);
            inputs.parse(obj,path_night,path_calibration,varargin{:});
            
            % reset
            obj.folder_data = [];
            obj.data_id      = [];
            obj.fitsHeader = [];
                        
            flagStatus = checkDataFolder(obj,path_night);
            if strcmp(flagStatus,'error')
                return
            else
                flagStatus = getDataID(obj);
            end
            
            if strcmp(flagStatus,'error')
                return
            else
                flagStatus = getImageHeader(obj);
            end
            
            if strcmp(flagStatus,'error')
                return
            else
                checkCalibrationFolder(obj,path_calibration);
            end
            obj.readingStatus = true;
        end
        
        %% Data loading
        function grabData(obj,varargin)
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'puako'));
            inputs.addParameter('id', 'ALL',@(x) ischar(x) | isstring(x) | iscell(x));
            inputs.addParameter('path_night', obj.folder_data.path_root,@ischar);
            inputs.addParameter('path_calibration',obj.folder_data.calibration, @ischar);
            inputs.parse(obj,varargin{:});
            
            %1\ Run the folder analysis if the user did not do it already
            id = inputs.Results.id;
            path_night =  inputs.Results.path_night;
            path_calibration = inputs.Results.path_calibration;
            if isempty(obj.folder_data.trs)
                readFolder(obj,path_night,path_calibration)
            end
            
            %2\ Select files the user want to process
            [id_trs,id_imag,fitsHdr,flagStatus] = selectFilesFromUser(obj,id);
            
            if strcmp(flagStatus,'error')
                fprintf('So I discontinue the data loading.\n Please check the format of the file name you''ve provided\n');
                return
            end
            
            % 3\ Check the size in Mb
            flagStatus = checkMemory(obj,id_trs,id_imag);
            % If the user decides the memory usage is too important, we quit
            if strcmp(flagStatus,'error')
                fprintf('So I discontinue the data loading.\n Please check data_id.trs.size to select  a sub-sample of files that I can handle within the available %.1f Gb memory ',dataSize/1e3);
                return
            end
            
            % 4\ Instantiating the trs class
            nObj = max(numel(id_trs),numel(id_imag));
            if nObj ==1
                path_trs = replace([obj.folder_data.trs,id_trs.id],'//','/');
                path_imag = replace([obj.folder_data.imag,id_imag.id],'//','/');
                obj.trs = telemetry({id_trs.name},obj.folder_data,path_trs,path_imag,fitsHdr);
            else
                hwait = waitbar(0,'I''m grabbing data, this may take some time ...');               
                tic;
                for kObj = 1:nObj
                    path_trs = replace([obj.folder_data.trs,id_trs(kObj).id],'//','/');
                    path_imag = replace([obj.folder_data.imag,id_imag(kObj).id],'//','/');
                    obj.trs(kObj) = telemetry({id_trs(kObj).name},obj.folder_data,path_trs,path_imag,fitsHdr(kObj));                    
                    waitbar(kObj/nObj);
                end
                close(hwait);
                fprintf('Done in %.2g s\n',toc());
            end
        end
        
        %% PSF reconstruction
        function getRecPSF(obj,varargin)
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'puako'));
            inputs.addParameter('id', 'ALL',@(x) iscell(x) | ischar(x) | isstring(x));
            inputs.addParameter('path_night', obj.folder_data.path_root,@ischar);
            inputs.addParameter('path_calibration',obj.folder_data.calibration, @ischar);
            inputs.parse(obj,varargin{:});
            
            %1\ Run the folder analysis if the user did not do it already
            id = inputs.Results.id;
            path_night =  inputs.Results.path_night;
            path_calibration = inputs.Results.path_calibration;
            if isempty(obj.folder_data.trs)
                if isempty(path_night) || isempty(path_calibration)
                    uiwait(errordlg('You must provide valid paths for the data and calibration folders !'));
                    return;
                else
                    readFolder(obj,path_night,path_calibration)
                end
            end
            
            %2\ Select files the user want to process
            [id_trs,id_imag,fitsHdr] = selectFilesFromUser(obj,id);
            
            %3\ Loop on files: so far puako is redoing the processing of
            %files it has already processed, will be fixed later.
            obj.res.psfr = [];
            obj.res.trs = [];
            for kObj = 1:numel(id_trs)
                % Get the telemetry
                path_trs = replace([obj.folder_data.trs,id_trs(kObj).id],'//','/');
                path_imag = replace([obj.folder_data.imag,id_imag(kObj).id],'//','/');
                obj.trs = telemetry({id_trs(kObj).name},obj.folder_data,path_trs,path_imag,fitsHdr(kObj));
                % Perform the PSFR
                obj.psfr = psfReconstruction(obj.trs);
                % Save results
                obj.res.psfr(kObj).image = obj.trs.cam.frame;
                obj.res.psfr(kObj).rec = tools.crop(obj.psfr.psf.rec,obj.trs.cam.resolution);
                obj.res.trs(kObj).seeing = obj.trs.res.seeing;
                obj.res.trs(kObj).noise = obj.trs.res.noise;
                obj.res.trs(kObj).zernike = obj.trs.res.zernike;
            end
        end
    end
end




