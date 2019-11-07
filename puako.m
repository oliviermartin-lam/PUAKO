classdef puako < handle
    
    properties
        % Data ID
        folder_data = [];                                         % Data path structure
        data_id = [];                                               % Data ID structure
        fitsHeader = [];                                           % Fits file header
        % Linked sub-classes
        trs;                                                             % telemetry reading and processing
        psfr;                                                            % Forward psf Reconstruction
        prime;                                                       % Hybrid psf reconstruction
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
            
            
            % Initialization
            obj.folder_data.trs = [];
            obj.folder_data.imag = [];
            obj.folder_data.date = [];
            obj.folder_data.path_root = [];
            obj.folder_data.calibration = [];
            obj.folder_data.massdimm = [];
            obj.folder_data.savings = [];
            
        end
        
        
        %% Data reading
        function readFolder(obj,path_night,path_calibration,varargin)
            inputs = inputParser;
            inputs.addRequired('obj', @(x) isa(x,'puako'));
            inputs.addRequired('path_night', @ischar);
            inputs.addRequired('path_calibration', @ischar);
            inputs.parse(obj,path_night,path_calibration,varargin{:});
            
            % Init
            obj.folder_data.trs = [];
            obj.folder_data.imag = [];
            obj.folder_data.date = [];
            obj.folder_data.path_root = [];
            obj.folder_data.calibration = [];
            obj.folder_data.massdimm = [];
            obj.folder_data.savings = [];
            
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
                obj.trs = [];
                for kObj = 1:nObj
                    path_trs = replace([obj.folder_data.trs,id_trs(kObj).id],'//','/');
                    path_imag = replace([obj.folder_data.imag,id_imag(kObj).id],'//','/');
                    if kObj == 1
                        obj.trs = telemetry({id_trs(kObj).name},obj.folder_data,path_trs,path_imag,fitsHdr(kObj));
                    else
                        obj.trs(kObj) = telemetry({id_trs(kObj).name},obj.folder_data,path_trs,path_imag,fitsHdr(kObj));
                    end
                    waitbar(kObj/nObj);
                end
                close(hwait);
                t0 = toc();
                
                if t0 < 60
                    fprintf('Done in %.2g s\n',t0);
                else
                    fprintf('Done in %.2g mn\n',t0/60);
                end
            end
        end
        
        %% PSF reconstruction
        function getRecPSF(obj,varargin)
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'puako'));
            inputs.addParameter('id', 'ALL',@(x) iscell(x) | ischar(x) | isstring(x));
            inputs.addParameter('path_night', obj.folder_data.path_root,@ischar);
            inputs.addParameter('path_calibration',obj.folder_data.calibration, @ischar);
            inputs.addParameter('path_save',obj.folder_data.savings, @ischar);
            inputs.parse(obj,varargin{:});
            
            %1\ Run the folder analysis if the user did not do it already
            id = inputs.Results.id;
            path_night =  inputs.Results.path_night;
            path_calibration = inputs.Results.path_calibration;
            path_save = inputs.Results.path_save;
            
            if isempty(obj.folder_data.trs) % data have not been loaded
                if isempty(path_night) || isempty(path_calibration)
                    uiwait(errordlg('You must provide valid paths for the data and calibration folders !'));
                    return;
                else
                    readFolder(obj,path_night,path_calibration);
                    
                    if ~isempty(path_save)
                        if ~isfolder(path_save)
                            fprintf(['I create the folder ',path_save,'t save results\n']);
                            mkdir(path_save);
                        end
                        obj.folder_data.savings = path_save;
                    end
                end
                
            else % some data have not been loaded already, but check whether the folder is the same
                tmp_fold = obj.folder_data;                
                flagStatus = checkDataFolder(obj,path_night);
                if strcmp(flagStatus,'error')
                    return;
                end
                % Reload the folder if necessary
                if isequal(obj.folder_data,tmp_fold)
                    fprintf(['I keep going with the same data folder.\n']);
                     if ~isequal(obj.data_id,obj.data_id)
                         fprintf(['But I redefine the data identities.\n']);
                         % Redefine the data identity
                         readFolder(obj,path_night,path_calibration)
                     end                    
                else
                    readFolder(obj,path_night,path_calibration);
                end
                
                % Create a savings folder
                if ~isempty(path_save)
                    if ~isfolder(path_save)
                        fprintf(['I create the folder ',path_save,'t save results\n']);
                        mkdir(path_save);
                    end
                    obj.folder_data.savings = path_save;
                end
                
            end
            
            tpsfr = tic();
            %2\ Select files the user want to process
            [id_trs,id_imag,fitsHdr] = selectFilesFromUser(obj,id);
            
            %3\ Loop on files: so far puako is redoing the processing of
            %files it has already processed, will be fixed later.
            nObj = numel(id_trs);
            obj.psfr = [];
            for kObj = 1:nObj
                % Get the telemetry
                path_trs = replace([obj.folder_data.trs,id_trs(kObj).id],'//','/');
                path_imag = replace([obj.folder_data.imag,id_imag(kObj).id],'//','/');
                obj.trs = telemetry({id_trs(kObj).name},obj.folder_data,path_trs,path_imag,fitsHdr(kObj));
                P             = obj.trs.tel.pupil;
                lambda   = obj.trs.cam.wavelength;
                Samp     = obj.trs.cam.samp;
                psInMas = obj.trs.cam.pixelScale;
                % Perform the PSFR
                if kObj == 1
                    obj.psfr = psfReconstruction(obj.trs);
                else
                    obj.psfr(kObj) = psfReconstruction(obj.trs);
                end
                %Save telemetry processing results
                hdr.r0 = obj.trs.res.seeing.r0;
                hdr.dr0 = obj.trs.res.seeing.dr0;
                hdr.L0 = obj.trs.res.seeing.L0;
                hdr.dL0 = obj.trs.res.seeing.dL0;
                hdr.w0 = obj.trs.res.seeing.w0;
                hdr.dw0 = obj.trs.res.seeing.dw0;
                hdr.seeing = obj.trs.res.seeing.seeing;
                hdr.dseeing = obj.trs.res.seeing.dseeing;
                hdr.mnoiseho = sqrt(obj.trs.res.noise.varn_ho)*1e9;
                hdr.mnoisett = sqrt(obj.trs.res.noise.varn_tt)*1e9;
                hdr.std_tip= std(obj.trs.tipTilt.slopes(1,:),[],2)*1e9;
                hdr.std_tilt= std(obj.trs.tipTilt.slopes(1,:),[],2)*1e9;
                hdr.std_ho      = sqrt(sum(std(obj.trs.rec.res(obj.trs.dm.pupilMask,:),[],2).^2)/nnz(obj.trs.dm.pupilMask))*1e9;
                
                % Get PSF results
                obj.psfr(kObj).sky = psfStats(obj.trs.cam.frame,P,lambda,Samp,psInMas,'flagMoffat',true,'flagGaussian',true);
                recpsf = tools.crop(obj.psfr(kObj).rec_,size(obj.trs.cam.frame));
                obj.psfr(kObj).psf  = psfStats(recpsf,P,lambda,Samp,psInMas,'im_ref',obj.trs.cam.frame);
                hdr.SRrec = obj.psfr(kObj).psf.SR;
                hdr.FWHMx = obj.psfr(kObj).psf.FWHMx;
                hdr.FWHMy = obj.psfr(kObj).psf.FWHMy;
                hdr.dFWHM = obj.psfr(kObj).psf.dFWHM;
                hdr.SRsky = obj.psfr(kObj).sky.SR;
                hdr.dSRsky = obj.psfr(kObj).sky.dSR;
                hdr.FWHMxsky = obj.psfr(kObj).sky.FWHMx;
                hdr.FWHMysky = obj.psfr(kObj).sky.FWHMy;
                hdr.dFWHMsky = obj.psfr(kObj).sky.dFWHM;
                
                % Get the error breakdown
                wfe = errorBreakDown( obj.psfr(kObj));
                
                hdr.sr_tot      = wfe.sr_tot;
                hdr.wfetot     = wfe.wfe_tot;
                hdr.wfencpa    = wfe.wfe_ncpa;
                hdr.wfefit     = wfe.wfe_fit;
                hdr.wfelag      = wfe.wfe_lag;
                hdr.wfenoise   = wfe.wfe_noise;
                hdr.wfealias   = wfe.wfe_alias;
                hdr.wfeaniso   = wfe.wfe_aniso;
                hdr.wfett  = wfe.wfe_tt;
                hdr.wfenoiTT = wfe.wfe_noiseTT;
                
                % Get precision
                hdr.gauss_x = obj.psfr(kObj).sky.catalogs.gaussian.x;
                hdr.gauss_y = obj.psfr(kObj).sky.catalogs.gaussian.y;
                hdr.gauss_dx = obj.psfr(kObj).sky.catalogs.gaussian.dx;
                hdr.gauss_dy = obj.psfr(kObj).sky.catalogs.gaussian.dy;
                hdr.gauss_f = obj.psfr(kObj).sky.catalogs.gaussian.flux;
                hdr.gauss_df = obj.psfr(kObj).sky.catalogs.gaussian.dflux;
                hdr.fvugauss = obj.psfr(kObj).sky.catalogs.gaussian.fvu;
                
                hdr.moff_x = obj.psfr(kObj).sky.catalogs.moffat.x;
                hdr.moff_y = obj.psfr(kObj).sky.catalogs.moffat.y;
                hdr.moff_dx = obj.psfr(kObj).sky.catalogs.moffat.dx;
                hdr.moff_dy = obj.psfr(kObj).sky.catalogs.moffat.dy;
                hdr.moff_f = obj.psfr(kObj).sky.catalogs.moffat.flux;
                hdr.moff_df = obj.psfr(kObj).sky.catalogs.moffat.dflux;
                hdr.fvumoff = obj.psfr(kObj).sky.catalogs.moffat.fvu;
                
                hdr.rec_x = obj.psfr(kObj).psf.catalogs.ref.x;
                hdr.rec_y = obj.psfr(kObj).psf.catalogs.ref.y;
                hdr.rec_dx = obj.psfr(kObj).psf.catalogs.ref.dx;
                hdr.rec_dy = obj.psfr(kObj).psf.catalogs.ref.dy;
                hdr.rec_f = obj.psfr(kObj).psf.catalogs.ref.flux;
                hdr.rec_df = obj.psfr(kObj).psf.catalogs.ref.dflux;
                hdr.fvurec = obj.psfr(kObj).psf.catalogs.ref.fvu;
                
                %Write results
                if ~isempty(obj.folder_data.savings)
                    filename = [obj.folder_data.savings,obj.trs.cam.name,'_',obj.trs.obj_name,'_psfr.fits'];
                    fits_write(filename,obj.psfr(kObj).psf.image,hdr);
                end
                
                if nObj > 1
                    % make sure that we do not store all the telemetry when looping over data
                    obj.psfr(kObj).trs.wfs.slopes = [];
                    obj.psfr(kObj).trs.dm.com = [];
                    obj.psfr(kObj).trs.rec = [];
                    obj.psfr(kObj).trs.mat = [];
                end
            end
            fprintf('Job done in %.2g s\n',toc(tpsfr));
        end
        
        %% Results Vizualization and statistical comparison
        
        function vizualizeResults(obj,varargin)
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'puako'));
            inputs.addParameter('resultsof',[], @ischar);
            inputs.addParameter('path_save',[], @ischar);
            inputs.parse(obj,varargin{:});
            
            resultsof = inputs.Results.resultsof;
            
            switch resultsof
                case 'psfr'
                    if isempty(obj.psfr)
                        fprintf('Sorry, you must run the PSF-R processing first, please use my getRecPSF facility\n');
                    else
                        displayResults(obj.psfr);
                    end
                case 'prime'
                    if isempty(obj.psfr)
                        fprintf('Sorry, you must run the PRIME processing first, please use my getPrimePSF facility\n');
                    else
                        displayResults(obj.prime);
                    end
                otherwise
                    fprintf('Sorry, you must provide ''PSF-R'' or ''PRIME'' for the ''resultsof'' entry\n');
            end
            
            
        end
    end
end




