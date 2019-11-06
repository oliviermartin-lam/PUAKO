classdef telemetry < handle
    
    properties (SetObservable=true)
        
        % ---------------------------- PATHS ---------------------------- %
        obj_name;
        date;
        path_trs;
        path_imag;
        path_calibration;
        path_massdimm;
        fitsHdr;
        % ------------------- SYSTEM/OBSERVING CONFIGURATION -------------------- %
        aoMode;
        sysTime;
        % Structures (not OOMAO class)
        src;                                                            % structure   containing science source info
        ngs;                                                            % structure  containing NGS info
        lgs;                                                              % structure  containing LGS info
        tel;                                                              % structure containing telescope info
        atm;                                                            % structure containing the atmosphere parameters
        wfs;                                                            %structure containing HO WFS info
        tipTilt;                                                        %structure containing TT info
        dm;                                                             %structure containing HO DM info
        cam;                                                           % structure containing Science camera info
        rec;                                                            % structure containing reconstructed wavefront
        mat;                                                           % structure containing system matrices
        holoop;                                                       % structure containing HO loop configuration
        ttloop;                                                        % structure containing TT loop configuration       
        res = [];                                                            % structure containing processing results
    end
    
    methods
        
        function obj = telemetry(obj_name,folder_data,path_trs,path_imag,fitsHdr)
            inputs = inputParser;           
            inputs.addRequired('obj_name', @iscell);
            inputs.addRequired('folder_data', @isstruct);
            inputs.addRequired('path_trs', @ischar);
            inputs.addRequired('path_imag', @ischar);
            inputs.addRequired('fitsHdr', @isstruct);
            inputs.parse(obj_name,folder_data,path_trs,path_imag,fitsHdr);
            
            %1\ Parsing inputs
            obj.obj_name = obj_name{1};
            obj.date         = inputs.Results.folder_data.date;
            obj.path_trs    = inputs.Results.path_trs;
            obj.path_imag      = inputs.Results.path_imag;
            obj.path_calibration   = inputs.Results.folder_data.calibration;
            obj.path_massdimm = inputs.Results.folder_data.massdimm;
            obj.fitsHdr = fitsHdr;
                  
            %2\ Initialize system structures
            obj = initializeStructures(obj);
            
            %3\ Restoring Calibrated data            
            obj = restoreCalibratedData(obj);
            
            %3\ Restoring AO telemetry            
            obj = restoreKeckTelemetry(obj);
                       
            %4\ Model temporal transfer functions
            obj = modelTransferFunctions(obj);
            
            %5\ Restoring and processing NIRC2 images
            obj = restoreDetectorImage(obj);
            
            %6\ Getting MASS/DIMM data         
            obj = restoreMassDimmMaunaKea(obj);           
            
        end      
        
        %% Zernike decomposition
        function std_n = getZernikeVariance(obj,jindex)
            if ~exist('jindex')
                jindex = 2:120;
            end
            inputs = inputParser;
            inputs.addRequired('obj', @(x) isa(x,'telemetry'));
            inputs.addRequired('jindex', @isnumeric);
            std_n= getZernikeDecomposition(obj,jindex);
        end
        
        %% Transfer function
        function displayAoTransferFunction(obj)           
            displayTransferFunction(obj);
        end
        
         %% Noise standard-deviation in nm
        function [stdn_ho,stdn_tt] = getNoiseSTD(obj,method)
            if nargin < 2
                method = 'autocorrelation';
            end
            if isstring(method) || ischar(method)
                method = {method};
            end
                     
            nObj = numel(obj);
            nMethod = numel(method);
            
            stdn_ho = zeros(nMethod,nObj);
            stdn_tt = zeros(nMethod,nObj);            
            for kObj=1:nObj
                obj(kObj).res.noise = [];
                for jM = 1:nMethod
                    tmp = estimateNoiseCovarianceFromTelemetry(obj(kObj),'method',method{jM});
                    stdn_ho(jM,kObj) = sqrt(tmp.varn_ho)*1e9;
                    stdn_tt(jM,kObj) = sqrt(tmp.varn_tt)*1e9;
                    obj(kObj).res.noise(jM).method = method{jM};
                    obj(kObj).res.noise(jM).varn_ho = stdn_ho(jM,kObj);
                    obj(kObj).res.noise(jM).varn_tt = stdn_tt(jM,kObj);
                    obj(kObj).res.noise(jM).Cn_ho   = tmp.Cn_ho;
                    obj(kObj).res.noise(jM).Cn_tt    = tmp.Cn_tt;
                end
            end
        end
        
         %% Seeing estimation
        function getSeeing(obj,varargin)
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'telemetry'));
            inputs.addParameter('nMin',4,@isnumeric);
            inputs.addParameter('nMax',120,@isnumeric);
            inputs.addParameter('fitL0',true,@islogical);
            inputs.addParameter('best',false,@islogical);
            inputs.addParameter('wvl',0.5e-6,@isnumeric);
            inputs.addParameter('aoMode','NGS',@ischar);
            inputs.addParameter('D1',9,@isnumeric);
            inputs.addParameter('D2',2.65,@isnumeric);
            inputs.parse(obj,varargin{:});
            

            nObj = numel(obj);                                              
            r0 = zeros(1,nObj);
            L0 = zeros(1,nObj);
            fwhm = zeros(1,nObj);
            dr0 = zeros(1,nObj);
            dL0 = zeros(1,nObj);
            dfwhm = zeros(1,nObj);
                                  
            for kObj=1:nObj
                    obj(kObj).res.seeing = [];
                    tmp = estimateSeeingFromTelemetry(obj(kObj),varargin{:});
                    r0(kObj) = tmp.r0;
                    L0(kObj) = tmp.L0;
                    fwhm(kObj) = tmp.fwhm;
                    dr0(kObj) = tmp.dr0;
                    dL0(kObj) = tmp.dL0;
                    dfwhm(kObj) = tmp.dfwhm;                    
                    obj(kObj).res.seeing = tmp;
            end                       
        end
        
    end
end

