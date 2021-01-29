classdef puako < handle
    
    properties
        % Data ID
        data_folder_path = [];                                % Data folder structure
        data_folder_id = [];                                    % Data ID structure
        % Linked sub-classes
        trs;                                                             % telemetry reading and processing
        psfr;                                                            % Forward psf Reconstruction
        psfp;                                                       % Hybrid psf reconstruction
        yesToAll;
        noToAll;
        threshSteeringMirror;
        keepInMemory = true;
    end
    
    properties (Dependent, SetObservable=true)
        path_imag;
        path_trs;
        path_calibration;
        file_fits;
        file_sav;
    end
    
    properties (Access=private)
        p_imag;
        p_trs;
        p_calib;
        p_fits;
        p_sav;
    end
    
    properties (Hidden = true)
        flagStatus;
        p_idx;
        trs_idx;
        imag_idx;
        file_idx;
    end
    
    methods
        
        function obj = puako(varargin)
            inputs = inputParser;
            inputs.addParameter('path_imag', [], @ischar);
            inputs.addParameter('path_trs', [], @ischar);
            inputs.addParameter('path_calibration', [], @ischar);
            inputs.addParameter('path_dark', [], @ischar);
            inputs.addParameter('path_sky', [], @ischar);
            inputs.addParameter('path_div', [], @ischar);
            inputs.addParameter('yesToAll', true, @islogical);
            inputs.addParameter('noToAll', false, @islogical);
            inputs.addParameter('threshSteeringMirror', 20, @isnumeric);
            inputs.parse(varargin{:})
            obj.yesToAll = inputs.Results.yesToAll;
            obj.noToAll = inputs.Results.noToAll;
            obj.threshSteeringMirror = inputs.Results.threshSteeringMirror;
            
            fprintf('\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
            fprintf('\t--\t\t Hello there, I''m PUAKO and I shall be your guide in this cruel AO world !\t\t\t--\n\n')
            fprintf(['\t--\t You can request me to do whatever is in my capabilities (lucky you), which are: \t\t\t--\n'...
                '\t--\t 1. Data reading\t: load IMAG/TRS files within a night folder and read fits header\t\t--\n'...
                '\t--\t 2. Data processing\t: process telemetry and image data\t\t\t\t--\n'...
                '\t--\t 3. PSF reconstruction\t: reconstruct the PSF in a forward mode (no best-fitting)\t\t\t--\n'...
                '\t--\t 4. Hybrid PSF-R\t: hybrid PSF reconstruction using pupil anf focal plane measurements\t\t--\n'...
                '\t--\t 5. Data analysis\t: display tools for analyzing AO performance and reconstruction results\t\t--'])
            fprintf('\n\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
            
            
            %Check folders
            if ~isempty(inputs.Results.path_dark)
                obj.data_folder_path.dark = inputs.Results.path_dark;
                if strcmp(obj.flagStatus,'error')
                    return
                end
            else
                obj.data_folder_path.dark = [];
            end
            
             if ~isempty(inputs.Results.path_sky)
                obj.data_folder_path.sky = inputs.Results.path_sky;
                if strcmp(obj.flagStatus,'error')
                    return
                end
             else
                 obj.data_folder_path.sky = [];
             end
            
             if ~isempty(inputs.Results.path_div)
                 obj.data_folder_path.div = inputs.Results.path_div;
                 if strcmp(obj.flagStatus,'error')
                     return
                 end
             else
                 obj.data_folder_path.div = [];
             end
              
            if ~isempty(inputs.Results.path_imag)
                obj.path_imag = inputs.Results.path_imag;
                if strcmp(obj.flagStatus,'error')
                    return
                end
            end
            
            
            if ~isempty(inputs.Results.path_trs)
                obj.path_trs = inputs.Results.path_trs;
                if strcmp(obj.flagStatus,'error')
                    return
                end
            end
            
            if ~isempty(inputs.Results.path_calibration)
                obj.path_calibration  = inputs.Results.path_calibration;
                if strcmp(obj.flagStatus,'error')
                    return
                end
            end
        end
        %% Set functions
        
        function set.path_imag(obj,val)
            % read the folder
            obj.flagStatus = readIMAGFolder(obj,val);
            if strcmp(obj.flagStatus,'error')
                return
            end
            obj.p_imag = val;
            % read the data
            obj.flagStatus = getIMAGDataID(obj);
            if strcmp(obj.flagStatus,'error')
                return
            end
        end
        function set.path_trs(obj,val)
            % read the folder
            obj.flagStatus = readTRSFolder(obj,val);
            if strcmp(obj.flagStatus,'error')
                return;
            end
            obj.p_trs = val;
            % read the data
            obj.flagStatus = getTRSDataID(obj);
            if strcmp(obj.flagStatus,'error')
                return
            end
        end
        function set.path_calibration(obj,val)
            obj.flagStatus = readCalibrationFolder(obj,val);
            if strcmp(obj.flagStatus,'error')
                return;
            end
            obj.p_calib = val;
        end
        function set.file_fits(obj,val)
            [obj.p_fits,obj.flagStatus ]= checkFitsFile(obj,val);
        end
        function set.file_sav(obj,val)
            [obj.p_sav,obj.flagStatus] = checkSavFile(obj,val);
        end
        
        %% Get functions
        function val = get.path_imag(obj)
            val = obj.p_imag;
        end
        function val = get.path_trs(obj)
            val = obj.p_trs;
        end
        function val = get.path_calibration(obj)
            val = obj.p_calib;
        end
        function val = get.file_fits(obj)
            val = obj.p_fits;
        end
        function val = get.file_sav(obj)
            val = obj.p_sav;
        end
        
        %% Data loading
        function grabData(obj,varargin)
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'puako'));
            inputs.addParameter('objname', {''},@(x) iscell(x) | ischar(x));
            inputs.addParameter('path_imag', obj.data_folder_path.imag, @ischar);
            inputs.addParameter('path_trs', obj.data_folder_path.trs, @ischar);
            inputs.addParameter('path_calibration',obj.data_folder_path.calibration, @ischar);
            inputs.addParameter('path_ncpa', '', @ischar);
            inputs.addParameter('path_save','', @ischar);
            inputs.addParameter('yesToAll', true, @islogical);
            inputs.addParameter('threshSteeringMirror', 20, @isnumeric);

            inputs.addParameter('resolution', 150, @isnumeric);
            inputs.addParameter('flagNoisemethod', 'autocorrelation', @ischar);
            inputs.addParameter('badModesList',[],@isnumeric);
            inputs.addParameter('jMin',4,@isnumeric);
            inputs.addParameter('jMax',120,@isnumeric);
            inputs.addParameter('fitL0',true,@islogical);
            inputs.addParameter('flagBest',false,@islogical);
            inputs.addParameter('flagMedian',false,@islogical);
            inputs.addParameter('D1',9,@isnumeric);
            inputs.addParameter('D2',2.65,@isnumeric);
            inputs.addParameter('getImageOnly',false,@islogical);
            inputs.addParameter('flagGaussian',true,@islogical);            
            inputs.addParameter('flagMoffat',true,@islogical);      
            inputs.parse(obj,varargin{:});

            resolution = inputs.Results.resolution;
            noisemethod = inputs.Results.flagNoisemethod;
            fitL0 = inputs.Results.fitL0;
            flagBest = inputs.Results.flagBest;
            flagMedian =  inputs.Results.flagMedian;
            D1 =inputs.Results.D1;
            D2 =inputs.Results.D2;
            badModesList =inputs.Results.badModesList;
            jMin = inputs.Results.jMin;
            jMax = inputs.Results.jMax;
            getImageOnly = inputs.Results.getImageOnly;
            flagMoffat = inputs.Results.flagMoffat;
            flagGaussian = inputs.Results.flagGaussian;
            path_ncpa = inputs.Results.path_ncpa;
            
            %1\ Select subsamples with memory check
            t0 = tic();                
            selectSamples(obj,inputs.Results.objname,'path_imag',inputs.Results.path_imag,'path_trs',inputs.Results.path_trs,'path_calibration',inputs.Results.path_calibration,...
                'path_save',inputs.Results.path_save,'yesToAll',inputs.Results.yesToAll,'threshSteeringMirror',inputs.Results.threshSteeringMirror,'getImageOnly',getImageOnly);
            if strcmp(obj.flagStatus,'error')
                return
            end
            
            % 2\ Instantiating the trs class
            obj.trs = [];
            if ~getImageOnly
                nObj = numel(obj.file_sav);
            else
                nObj = numel(obj.file_fits);
            end
            
            fprintf('I''m grabbing data, this may take some time ...');
            for kObj = 1:nObj
                %2.1 Get the appropriate TRS/IMAG paths and header
                pathIMAG = obj.file_fits{kObj};
                if ~getImageOnly
                    pathTRS = obj.file_sav{kObj};
                    hdr_k = obj.data_folder_id.trs(obj.trs_idx(kObj)).hdr;
                    objname = {obj.data_folder_id.trs(obj.trs_idx(kObj)).name};
                else
                    pathTRS = '';
                    hdr_k = obj.data_folder_id.imag(obj.imag_idx(kObj)).hdr;
                    objname = {obj.data_folder_id.imag(obj.imag_idx(kObj)).name};
                end
                %2.2 Get the telemetry
                if nObj == 1
                    obj.trs = telemetry(objname,pathTRS,pathIMAG,obj.data_folder_path,hdr_k,'path_ncpa',path_ncpa,'resolution',resolution,'flagNoisemethod',noisemethod,'fitL0',fitL0,'flagBest',flagBest,'flagMedian',...
                flagMedian,'D1',D1,'D2',D2,'badModesList',badModesList,'jMin',jMin,'jMax',jMax,'getImageOnly',getImageOnly,'flagMoffat',flagMoffat,'flagGaussian',flagGaussian);
                else
                    obj.trs{kObj} = telemetry(objname,pathTRS,pathIMAG,obj.data_folder_path,hdr_k,'path_ncpa',path_ncpa,'resolution',resolution,'flagNoisemethod',noisemethod,'fitL0',fitL0,'flagBest',flagBest,'flagMedian',...
                flagMedian,'D1',D1,'D2',D2,'badModesList',badModesList,'jMin',jMin,'jMax',jMax,'getImageOnly',getImageOnly,'flagMoffat',flagMoffat,'flagGaussian',flagGaussian);
                end
            end            
            
            
            t0 = toc(t0);
            if t0 < 60
                fprintf('Done in %.2g s\n',t0);
            else
                fprintf('Done in %.2g mn\n',t0/60);
            end
        end
        
        %% PSF reconstruction
        %z = zernike_puako(4:5,p.trs.tel.resolution);p.psfr = psfReconstruction(p.trs,'fov',300,'statModes',{[300,300],z.modes});
        function getRecPSF(obj,varargin)
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'puako'));
            inputs.addParameter('objname',{''}, @(x) iscell(x) | ischar(x));
            inputs.addParameter('path_imag', obj.data_folder_path.imag, @ischar);
            inputs.addParameter('path_trs', obj.data_folder_path.trs, @ischar);
            inputs.addParameter('path_calibration',obj.data_folder_path.calibration, @ischar);
            inputs.addParameter('path_ncpa', '', @ischar);
            inputs.addParameter('path_save',[], @ischar);
            inputs.addParameter('yesToAll', true, @islogical);
            inputs.addParameter('keepInMemory', false, @islogical);
            inputs.addParameter('threshSteeringMirror', 20, @isnumeric);
            inputs.addParameter('resolution', 150, @isnumeric);
            inputs.addParameter('statModes', {[]}, @iscell);
            inputs.addParameter('flagAoPattern','circle', @ischar);            
            inputs.addParameter('flagDphimethod', 'Vii', @ischar);
            inputs.addParameter('fitBg',false, @islogical);
            inputs.addParameter('flagNoisemethod', 'autocorrelation', @ischar);
            inputs.addParameter('badModesList',[],@isnumeric);
            inputs.addParameter('jMin',4,@isnumeric);
            inputs.addParameter('jMax',120,@isnumeric);
            inputs.addParameter('fitL0',true,@islogical);
            inputs.addParameter('flagBest',false,@islogical);
            inputs.addParameter('flagMedian',false,@islogical);
            inputs.addParameter('D1',9,@isnumeric);
            inputs.addParameter('D2',2.65,@isnumeric);
            inputs.addParameter('fov',[],@isnumeric);
            inputs.parse(obj,varargin{:});
            
            %1\ Checking if the user wants to save results or not
            obj.data_folder_path.savings = inputs.Results.path_save;
            if isempty(obj.data_folder_path.savings) && ~obj.yesToAll
                msg = input('Your results are not going to be saved, are you sure you want to continue (''y''/''n'') ?');
                if contains(upper(msg),'N')
                    fprintf('Ok I quit, please redo the call and provide an existing saving path by using the path_save property\n');
                    return
                end
                fprintf('Ok here we go.\n');
            end
            
            %2\ Getting user's inputs
            pathimg             = inputs.Results.path_imag;
            pathtrs             = inputs.Results.path_trs;
            pathcalibration     = inputs.Results.path_calibration;
            path_ncpa           = inputs.Results.path_ncpa;
            yToAll              = inputs.Results.yesToAll;
            obj.keepInMemory    = inputs.Results.keepInMemory;
            threshSteer         = inputs.Results.threshSteeringMirror;
            resolution          = inputs.Results.resolution;
            statModes           = inputs.Results.statModes;
            flagAoPattern       = inputs.Results.flagAoPattern;
            flagDphimethod      = inputs.Results.flagDphimethod;
            fitBg               = inputs.Results.fitBg;
            noisemethod         = inputs.Results.flagNoisemethod;
            fitL0               = inputs.Results.fitL0;
            flagBest            = inputs.Results.flagBest;
            flagMedian          = inputs.Results.flagMedian;
            D1                  = inputs.Results.D1;
            D2                  = inputs.Results.D2;
            badModesList        = inputs.Results.badModesList;
            jMin                = inputs.Results.jMin;
            jMax                = inputs.Results.jMax;
            fov_fit             = inputs.Results.fov;
            if isempty(fov_fit)
                fov_fit         = 4/3*inputs.Results.resolution;
            end
            fov_fit             = max(resolution + 4,fov_fit);

            %3\ Select subsamples with memory check
            t0 = tic();
            obj.keepInMemory = inputs.Results.keepInMemory;
            selectSamples(obj,inputs.Results.objname,'path_imag',pathimg,'path_trs',...
                pathtrs,'path_calibration',pathcalibration,'yesToAll',yToAll,'threshSteeringMirror',threshSteer);
            if strcmp(obj.flagStatus,'error')
                return
            end
            
            %4\ Loop on corresponding files.
            obj.psfr = [];
            nObj = numel(obj.file_sav);
            fprintf(['I''m performing PSFR over ',num2str(nObj),' files, this may take some time ...']);
            for kObj = 1:nObj
                %4.1 Get the appropriate TRS/IMAG paths and header
                pathTRS     = obj.file_sav{kObj};
                pathIMAG    = obj.file_fits{kObj};
                hdr_k       = obj.data_folder_id.trs(obj.trs_idx(kObj)).hdr;
                objname     = {obj.data_folder_id.trs(obj.trs_idx(kObj)).name};
                %4.2 Grab data
                obj.trs = telemetry(objname,pathTRS,pathIMAG,obj.data_folder_path,hdr_k,'path_ncpa',path_ncpa,'resolution',resolution,'flagNoisemethod',noisemethod,'fitL0',fitL0,'flagBest',flagBest,'flagMedian',...
                flagMedian,'D1',D1,'D2',D2,'badModesList',badModesList,'jMin',jMin,'jMax',jMax);

                %4.3 PSF reconstruction
                if obj.keepInMemory % we increment obj.psfr
                    obj.psfr{kObj} = psfReconstruction(obj.trs,'statModes',statModes,'flagAoPattern',flagAoPattern,'fitBg',fitBg,'fov',fov_fit,'flagDphimethod',flagDphimethod);
                    im  = obj.psfr{kObj}.psf.image;
                else % we do not not increment obj.psfr so we save memory
                    obj.psfr = psfReconstruction(obj.trs,'statModes',statModes,'flagAoPattern',flagAoPattern,'fitBg',fitBg,'fov',fov_fit,'flagDphimethod',flagDphimethod);
                    im  = obj.psfr.psf.image;
                end
                
                %5\ Write results
                if ~isempty(obj.data_folder_path.savings)
                    hdr = defineHeaderFromResults(obj.psfr);
                    filename = [obj.data_folder_path.savings,obj.trs.cam.name,'_',obj.trs.obj_name,'_psfr.fits'];
                    fits_write(filename,im,hdr);
                end
            end
            
            %6\ Get the time spent
            t0 = toc(t0);
            if t0 < 60
                fprintf('Done in %.2g s\n',t0);
            else
                fprintf('Done in %.2g mn\n',t0/60);
            end
            
        end
        
        %% PRIME reconstruction
        function getPrimePSF(obj,varargin)
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'puako'));
            inputs.addParameter('objname', {''}, @(x) iscell(x) | ischar(x));
            inputs.addParameter('path_imag', obj.data_folder_path.imag, @ischar);
            inputs.addParameter('path_trs', obj.data_folder_path.trs, @ischar);
            inputs.addParameter('path_calibration',obj.data_folder_path.calibration, @ischar);
            inputs.addParameter('path_ncpa', '', @ischar);
            inputs.addParameter('path_save',[], @ischar);
            inputs.addParameter('yesToAll', false, @islogical);
            inputs.addParameter('keepInMemory', false, @islogical);
            inputs.addParameter('threshSteeringMirror', 20, @isnumeric);
            inputs.addParameter('flagBrightestStar', false, @islogical);
            inputs.addParameter('resolution', 150, @isnumeric);
            inputs.addParameter('fitR0',true,@islogical);
            inputs.addParameter('fitCn2',false,@islogical);
            inputs.addParameter('fitGains',[true,true,true],@islogical);
            inputs.addParameter('jZernGain',[],@isnumeric);
            inputs.addParameter('fitStatModes',[],@isnumeric);
            inputs.addParameter('x_fixed',{[]},@iscell);
            inputs.addParameter('fitBg',true,@islogical);
            inputs.addParameter('ron',0,@isnumeric);
            inputs.addParameter('flagJacobian',false,@islogical);
            inputs.addParameter('algorithm','Trust-region-reflective',@ischar);
            inputs.addParameter('statModesFunction','', @ischar);
            inputs.addParameter('flagAoPattern','circle', @ischar);            
            inputs.addParameter('flagNoisemethod', 'autocorrelation', @ischar);
            inputs.addParameter('flagDphimethod', 'Vii', @ischar);
            inputs.addParameter('badModesList',[],@isnumeric);
            inputs.addParameter('jMin',4,@isnumeric);
            inputs.addParameter('jMax',120,@isnumeric);
            inputs.addParameter('fitL0',true,@islogical);
            inputs.addParameter('flagBest',false,@islogical);
            inputs.addParameter('flagMedian',false,@islogical);
            inputs.addParameter('D1',9,@isnumeric);
            inputs.addParameter('D2',2.65,@isnumeric);
            inputs.addParameter('fov',[],@isnumeric);
            inputs.addParameter('umax',10,@isnumeric);
            inputs.addParameter('aoinit',10,@isnumeric);
            inputs.addParameter('aobounds',10,@isnumeric);
            inputs.parse(obj,varargin{:});
            
           %1\ Checking if the user wants to save results or not
            obj.data_folder_path.savings = inputs.Results.path_save;
            if isempty(obj.data_folder_path.savings) && ~obj.yesToAll
                msg = input('Your results are not going to be saved, are you sure you want to continue (''y''/''n'') ?');
                if contains(upper(msg),'N')
                    fprintf('Ok I quit, please redo the call and provide an existing saving path by using the path_save property\n');
                    return
                end
                fprintf('Ok here we go.\n');
            end
            
            %2\ Getting user's inputs
            
            % Data inputs
            pathimg         = inputs.Results.path_imag;
            pathtrs         = inputs.Results.path_trs;
            pathcalibration = inputs.Results.path_calibration;
            path_ncpa       = inputs.Results.path_ncpa;
            yToAll          = inputs.Results.yesToAll;
            obj.keepInMemory= inputs.Results.keepInMemory;
            threshSteer     = inputs.Results.threshSteeringMirror;         
            % Fitting options
            flagBrightestStar= inputs.Results.flagBrightestStar;
            flagJacobian    = inputs.Results.flagJacobian;
            jZernGain       = inputs.Results.jZernGain;
            fitR0           = inputs.Results.fitR0;
            fitGains        = inputs.Results.fitGains;
            fitCn2          = inputs.Results.fitCn2;
            fitBg           = inputs.Results.fitBg;
            ron             = inputs.Results.ron;
            fitStatModes    = inputs.Results.fitStatModes;
            x_fixed         = inputs.Results.x_fixed;            
            algo            = inputs.Results.algorithm;
            % PSFR options
            statModesFunction = inputs.Results.statModesFunction;            
            flagAoPattern   = inputs.Results.flagAoPattern;
            % Telemetry processing options
            noisemethod     = inputs.Results.flagNoisemethod;
            flagDphimethod  = inputs.Results.flagDphimethod;
            fitL0           = inputs.Results.fitL0;
            flagBest        = inputs.Results.flagBest;
            flagMedian      = inputs.Results.flagMedian;
            D1              = inputs.Results.D1;
            D2              = inputs.Results.D2;
            badModesList    = inputs.Results.badModesList;
            jMin            = inputs.Results.jMin;
            jMax            = inputs.Results.jMax;
            resolution      = inputs.Results.resolution;
            fov_fit         = inputs.Results.fov;
            umax            = inputs.Results.umax;
            aoinit          = inputs.Results.aoinit;
            aobounds        = inputs.Results.aobounds;
            if isempty(fov_fit)
                fov_fit     = 4/3*inputs.Results.resolution;
            end
            fov_fit         = max(resolution + 4,fov_fit);
            
            %3\ Select subsamples regarding user's input with memory check
            t0 = tic();
            obj.keepInMemory = inputs.Results.keepInMemory;
            selectSamples(obj,inputs.Results.objname,'path_imag',pathimg,'path_trs',...
                pathtrs,'path_calibration',pathcalibration,'yesToAll',yToAll,'threshSteeringMirror',threshSteer);
            if strcmp(obj.flagStatus,'error')
                return
            end
            
            %4\ Loop on corresponding files
            obj.psfp = [];
            nObj = numel(obj.file_sav);
            fprintf(['I''m performing PSFR and PRIME over ',num2str(nObj),' files, this may take some time ...']);
            for kObj = 1:nObj
                %4.1 Get the appropriate TRS/IMAG paths and header
                pathTRS     = obj.file_sav{kObj};
                pathIMAG    = obj.file_fits{kObj};
                hdr_k       = obj.data_folder_id.trs(obj.trs_idx(kObj)).hdr;
                objname     = {obj.data_folder_id.trs(obj.trs_idx(kObj)).name};
                
                %4.2 Get the telemetry
                obj.trs = telemetry(objname,pathTRS,pathIMAG,obj.data_folder_path,hdr_k,'path_ncpa',path_ncpa,'resolution',resolution,'flagNoisemethod',noisemethod,'fitL0',fitL0,'flagBest',flagBest,'flagMedian',...
                flagMedian,'D1',D1,'D2',D2,'badModesList',badModesList,'jMin',jMin,'jMax',jMax,'flagBrightestStar',flagBrightestStar,'umax',umax);

                %4.3 PSF reconstruction
                obj.psfr = psfReconstruction(obj.trs,'flagAoPattern',flagAoPattern,'fitBg',fitBg,'fov',fov_fit,'flagDphimethod',flagDphimethod);
                
                %4.4 Perform the Hybrid PSFR
                if obj.keepInMemory % we increment obj.psfp
                    obj.psfp{kObj} = prime(obj.psfr,'jZernGain',jZernGain,'ron',ron,...
                        'flagJacobian',flagJacobian,'fitGains',fitGains,'fitCn2',fitCn2,'x_fixed',x_fixed);
                    im  = obj.psfp{kObj}.psf.image;
                else % we don't increment obj.psfp so we save memory
                    obj.psfp = prime(obj.psfr,'algorithm',algo,'x_fixed',x_fixed,'ron',ron,'flagJacobian',flagJacobian,...
                        'fitR0',fitR0,'fitCn2',fitCn2,'fitGains',fitGains,'jZernGain',jZernGain,'fitBg',fitBg,...
                        'fitStatModes',fitStatModes,'statModesFunction',statModesFunction,'aoinit',aoinit,'aobounds',aobounds);
                    im  = obj.psfp.psf.image;                   
                end

                %5\ Write results
                if ~isempty(obj.data_folder_path.savings)
                    hdr = defineHeaderFromResults(obj.psfp);
                    filename = [obj.data_folder_path.savings,obj.trs.cam.name,'_',obj.trs.obj_name,'_prime.fits'];
                    fits_write(filename,im,hdr);
                end
            end
            
            %6\ Get the time spent
            t0 = toc(t0);
            if t0 < 60
                fprintf('Done in %.2g s\n',t0);
            else
                fprintf('Done in %.2g mn\n',t0/60);
            end
        end
        
        %% Results Vizualization and statistical comparison        
        function vizualizeResults(obj,varargin)
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'puako'));
            inputs.addParameter('resultsof',[], @ischar);
            inputs.addParameter('path_save',[], @ischar);
            inputs.addParameter('kObj',1, @isnumeric);
            inputs.addParameter('fov',obj.trs.cam.resolution, @isnumeric);
            inputs.parse(obj,varargin{:});
            
            resultsof = inputs.Results.resultsof;
            kObj = inputs.Results.kObj;
            fov = inputs.Results.fov;
            
            switch resultsof
                case 'psfr'
                    if isempty(obj.psfr)
                        fprintf('Sorry, you must run the PSF-R processing first, please use my getRecPSF facility\n');
                    else
                        if kObj > numel(obj.psfr)
                            fprintf('Sorry, I have''nt processed as many files.\n');
                        else
                            if numel(obj.psfr)>1
                                displayResults(obj.psfr{kObj},'fov',fov);
                            else
                                displayResults(obj.psfr,'fov',fov);
                            end
                        end
                    end
                case 'prime'
                    if isempty(obj.psfp)
                        fprintf('Sorry, you must run the PRIME processing first, please use my getPrimePSF facility\n');
                    else
                        if kObj > numel(obj.psfp)
                            fprintf('Sorry, I have''nt processed as many files.\n');
                        else
                            if numel(obj.psfp)>1
                                displayResults(obj.psfp{kObj},'fov',fov);
                            else
                                displayResults(obj.psfp,'fov',fov);
                            end
                        end
                    end
                otherwise
                    fprintf('Sorry, you must provide ''PSF-R'' or ''PRIME'' for the ''resultsof'' entry\n');
            end                        
        end
                         
    end
end




