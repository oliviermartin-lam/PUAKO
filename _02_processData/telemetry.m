classdef telemetry < handle
    
    properties (SetObservable=true)
        
        % ---------------------------- PATHS ---------------------------- %
        obj_name;
        date;
        path_trs;
        path_imag;
        path_calibration;
        path_ncpa;
        path_massdimm;
        path_dark;
        path_sky;
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
        res = [];                                                      % structure containing processing results
        sky;                                                            % psfStats class for the sky image
    end
    
    methods
        
        function obj = telemetry(obj_name,path_trs,path_imag,folder_paths,fitsHdr,varargin)
            inputs = inputParser;           
            inputs.addRequired('obj_name', @(x) isa(x,'cell') || isa(x,'aoSystem'));
            inputs.addRequired('path_trs', @ischar);
            inputs.addRequired('path_imag', @ischar);
            inputs.addRequired('folder_paths', @isstruct);
            inputs.addRequired('fitsHdr', @iscell);
            inputs.addParameter('path_ncpa', [], @ischar);
            inputs.addParameter('resolution',150,@isnumeric);
            inputs.addParameter('flagNoisemethod','autocorrelation',@ischar);
            inputs.addParameter('badModesList',[],@isnumeric);
            inputs.addParameter('jMin',4,@isnumeric);
            inputs.addParameter('jMax',120,@isnumeric);
            inputs.addParameter('fitL0',true,@islogical);
            inputs.addParameter('flagBest',false,@islogical);
            inputs.addParameter('flagMedian',false,@islogical);
            inputs.addParameter('wvl',0.5e-6,@isnumeric);
            inputs.addParameter('D1',9,@isnumeric);
            inputs.addParameter('D2',2.65,@isnumeric);            
            inputs.addParameter('getImageOnly',false,@islogical);            
            inputs.addParameter('flagGaussian',false,@islogical);            
            inputs.addParameter('flagMoffat',false,@islogical);      
            inputs.addParameter('flagBrightestStar',false,@islogical);
            inputs.addParameter('umax',10,@isnumeric);
            
            inputs.parse(obj_name,path_trs,path_imag,folder_paths,fitsHdr,varargin{:});
                                
            %1\ Initialization
            obj = initializeStructures(obj);
            obj.cam.resolution = inputs.Results.resolution;
            
            %2\ Parsing inputs
            fitL0           = inputs.Results.fitL0;
            flagBest        = inputs.Results.flagBest;
            flagMedian      = inputs.Results.flagMedian;
            D1              = inputs.Results.D1;
            D2              = inputs.Results.D2;
            badModesList    = inputs.Results.badModesList;
            jMin            = inputs.Results.jMin;
            jMax            = inputs.Results.jMax;
            noiseMethod     = inputs.Results.flagNoisemethod;
            getImageOnly    = inputs.Results.getImageOnly;
            flagMoffat      = inputs.Results.flagMoffat;
            flagGaussian    = inputs.Results.flagGaussian;
            flagBrightestStar = inputs.Results.flagBrightestStar;
            umax            = inputs.Results.umax;
            
            if isa(obj_name,'aoSystem')    
                % Interface with KASP
                obj     = fromAoSystemClassToTelemetry(obj,obj_name);      
                aoMode  = obj.aoMode;
            else
                obj.obj_name        = obj_name{1};
                obj.path_trs        = inputs.Results.path_trs;
                obj.path_imag       = inputs.Results.path_imag;
                obj.path_calibration= inputs.Results.folder_paths.calibration;
                obj.path_ncpa       = inputs.Results.path_ncpa;
                obj.path_massdimm   = inputs.Results.folder_paths.massdimm;
                obj.path_dark       = inputs.Results.folder_paths.dark;
                obj.path_sky        = inputs.Results.folder_paths.sky;
                obj.fitsHdr         = fitsHdr;                                              
                obj.date            = cell2mat(obj.fitsHdr(strcmp(obj.fitsHdr(:,1),'DATE-OBS'),2));
                obj.date            = strjoin(split(obj.date,'-'),'');
                
                % detect if the laser was turned on
                if strcmp(cell2mat(obj.fitsHdr(contains(obj.fitsHdr(:,1),'LSPROP'),2)),'no')
                    obj.aoMode = 'NGS';
                else
                    obj.aoMode = 'LGS';
                end
                    
                
                %3\ Restoring Calibrated data
                obj = restoreCalibratedData(obj,getImageOnly);
                
                %4\ Restoring and processing NIRC2 images
                obj = processDetectorImage(obj,'flagMoffat',flagMoffat,'flagGaussian',flagGaussian,'flagBrightestStar',flagBrightestStar,'umax',umax);
                
                if ~getImageOnly
                    %5\ Restoring AO telemetry
                    obj = restoreKeckTelemetry(obj);
                    
                    %6\ Getting MASS/DIMM data
                    obj = restoreMassDimmMaunaKea(obj);
                end
            end
            
            if ~getImageOnly
                %7\ Estimate the number of photons
                obj = estimateNumberPhotons(obj);
                
                %8\ Model temporal transfer functions
                obj = modelTransferFunctions(obj);
                
                %9\ Data processing: noise estimation
                obj.res.noise = estimateNoiseCovarianceFromTelemetry(obj,'flagNoisemethod',noiseMethod);
                
                %10\ Data processing: seeing estimation
                [obj.res.seeing, obj.res.zernike] = estimateSeeingFromTelemetry(obj,'fitL0',fitL0,'flagBest',flagBest,'flagMedian',...
                    flagMedian,'D1',D1,'D2',D2,'badModesList',badModesList,'aoMode',obj.aoMode,'jMin',jMin,'jMax',jMax);
            end
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
        
        function displayZernikeFitting(obj)
            inputs = inputParser;
            inputs.addRequired('obj', @(x) isa(x,'telemetry'));
                        
            if (isfield(obj,'res') || isprop(obj,'res')) && isfield(obj.res,'zernike') && ~isempty(obj.res.zernike)
                h = figure;
                semilogy(obj.res.zernike.jindex,obj.res.zernike.std_meas,'ks--','MarkerFaceColor','k','MarkerSize',5);
                hold on;
                semilogy(obj.res.zernike.jindex,obj.res.zernike.std_model ,'ro--','MarkerFaceColor','r','MarkerSize',5);
                semilogy(obj.res.zernike.jindex,obj.res.zernike.std_noise,'bd--','MarkerFaceColor','b','MarkerSize',5);
                xlabel('Noll''s j-index','interpreter','latex','FontSize',20);
                ylabel('Zernike coefficients std (nm)','interpreter','latex','FontSize',20);
                set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex');
                legend({'Measurements','Model-fitted','Noise'},'interpreter','latex','FontSize',20);
            else
                msg = input('Sorry, you must active the seeing estimation procedure first, do you want me to do it (''y/n'') ?');
                if strcmpi(msg,'Y')
                    obj.getSeeing();
                    displayZernikeFitting(obj);
                end
            end
        end
        
        %% Transfer function
        function displayAoTransferFunction(obj)           
            displayTransferFunction(obj);
        end
        
         %% Noise standard-deviation in nm
        function [stdn_ho,stdn_tt] = getNoiseSTD(obj,flagNoisemethod)
            if nargin < 2
                flagNoisemethod = 'autocorrelation';
            end
            if isstring(flagNoisemethod) || ischar(flagNoisemethod)
                flagNoisemethod = {flagNoisemethod};
            end
                     
            nObj = numel(obj);
            nMethod = numel(flagNoisemethod);
            
            stdn_ho = zeros(nMethod,nObj);
            stdn_tt = zeros(nMethod,nObj);            
            for kObj=1:nObj
                obj(kObj).res.noise = [];
                for jM = 1:nMethod
                    tmp = estimateNoiseCovarianceFromTelemetry(obj(kObj),'method',flagNoisemethod{jM});
                    stdn_ho(jM,kObj) = sqrt(tmp.varn_ho)*1e9;
                    stdn_tt(jM,kObj) = sqrt(tmp.varn_tt)*1e9;
                    obj(kObj).res.noise(jM).method = flagNoisemethod{jM};
                    obj(kObj).res.noise(jM).varn_ho = stdn_ho(jM,kObj);
                    obj(kObj).res.noise(jM).varn_tt = stdn_tt(jM,kObj);
                    obj(kObj).res.noise(jM).Cn_ho   = tmp.Cn_ho;
                    obj(kObj).res.noise(jM).Cn_tt    = tmp.Cn_tt;
                end
            end
        end
        
         %% Seeing estimation
        function out = getSeeing(obj,varargin)
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'telemetry'));
            inputs.addParameter('nMin',4,@isnumeric);
            inputs.addParameter('nMax',120,@isnumeric);
            inputs.addParameter('fitL0',true,@islogical);
            inputs.addParameter('flagBest',false,@islogical);
            inputs.addParameter('flagMedian',true,@islogical);
            inputs.addParameter('wvl',0.5e-6,@isnumeric);
            inputs.addParameter('aoMode','NGS',@ischar);
            inputs.addParameter('D1',11.25,@isnumeric);
            inputs.addParameter('D2',2.65,@isnumeric);
            inputs.parse(obj,varargin{:});
            

            nObj = numel(obj);                                              
            r0 = zeros(1,nObj);
            L0 = zeros(1,nObj);
            seeing = zeros(1,nObj);
            dr0 = zeros(1,nObj);
            dL0 = zeros(1,nObj);
            dseeing = zeros(1,nObj);
                                  
            for kObj=1:nObj
                    obj(kObj).res.seeing = [];
                    [out(kObj),resZ] = estimateSeeingFromTelemetry(obj(kObj),varargin{:});
                    r0(kObj) = out.r0;
                    L0(kObj) = out.L0;
                    seeing(kObj) = out.seeing;
                    dr0(kObj) = out.dr0;
                    dL0(kObj) = out.dL0;
                    dseeing(kObj) = out.dseeing;                        
                    obj(kObj).res.seeing = out;
                    obj(kObj).res.zernike = resZ;
            end                       
        end
        
    end
end

