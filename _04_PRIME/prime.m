classdef prime < handle
    
    properties
        % PSF Model
        psfr;
        statModes;
        % --------------------- Fitting process setup        
        fov_sky;
        fov_fit;
        ron;
        modelFUN;
        x_init;
        lbounds;
        ubounds;
        fitOption;
        ydata;
        xdata;
        weightMap;
        normFactor;
        jZernGain;                
        nGainsHO;
        idxR0;
        idxCn2;
        idxDao;
        idxDtt;
        idxDal;
        idxStatModes;
        fitBg;
        x_final;
        x_prec;
        x_fixed;
        list_fixed;
        beta_;
        JacobianMatrix_;
        residualMap_;
        rec_;
        % Stars Initial guess
        idSrc;
        nParamStars;            % Number of stars parameters
        initStars;              % Stars parameters initial guess
        xStars;yStars;          % Astrometry initial guess
        fluxStars;              % Photometry initial guess
        phasor;
        % --------------------- PSF outputs
        psf;
        % --------------------- Retrieved outputs
        atm_fit;
        gains_fit;
        catalog_fit;              
        map_fit;   
        flagError;
    end
    
    methods
        function obj = prime(psfr,varargin)
            inputs = inputParser;
            inputs.addRequired('psfr',@(x) isa(x,'psfReconstruction'));
            inputs.addParameter('aoinit',[],@isnumeric);
            inputs.addParameter('aobounds',[],@isnumeric);
            inputs.addParameter('MaxIter',500,@isnumeric);
            inputs.addParameter('TolX',1e-15,@isnumeric);
            inputs.addParameter('TolFun',1e-15,@isnumeric);
            inputs.addParameter('MaxFunEvals',5e3,@isnumeric);
            inputs.addParameter('InitDamping',1e-5,@isnumeric);
            inputs.addParameter('display','iter',@ischar);
            inputs.addParameter('ron',0,@isnumeric);
            inputs.addParameter('x_fixed',{[]},@iscell);
            inputs.addParameter('fitR0',true,@islogical);
            inputs.addParameter('fitCn2',false,@islogical);
            inputs.addParameter('fitGains',[true,true,false],@islogical);
            inputs.addParameter('fitBg',true,@islogical);
            inputs.addParameter('statModesFunction','',@ischar);
            inputs.addParameter('fitStatModes',[],@isnumeric);
            inputs.addParameter('jZernGain',[],@isnumeric);
            inputs.addParameter('flagJacobian',false,@islogical);
            inputs.addParameter('algorithm','trust-region-reflective ',@ischar);
            inputs.addParameter('idSrc',1,@isnumeric);
            inputs.addParameter('umax',5,@isnumeric);
            inputs.parse(psfr,varargin{:});

            % Parse inputs            
            obj.psfr     = psfr;
            MaxIter      = inputs.Results.MaxIter;
            TolX         = inputs.Results.TolX;
            TolFun       = inputs.Results.TolFun;
            algo         = inputs.Results.algorithm;
            MaxFunEvals  = inputs.Results.MaxFunEvals;
            InitDamping  = inputs.Results.InitDamping;
            display      = inputs.Results.display;
            obj.ron      = inputs.Results.ron;
            flagJacobian = inputs.Results.flagJacobian;
            obj.list_fixed = inputs.Results.x_fixed;
            obj.fitBg    = inputs.Results.fitBg;
            obj.idSrc    = inputs.Results.idSrc;
            
            %1\ Check the image dimensions
            obj.fov_sky  = size(obj.psfr.trs.cam.image,1);
            obj.fov_fit  = size(obj.psfr.rec_,1);
            
            %2\ Normalize the observation and define the weight matrix
            obj.ydata       = obj.psfr.trs.cam.image(:,:,obj.idSrc);
            obj.normFactor  = sum(obj.ydata(obj.ydata>0));

            %3\ Define the weighting matrix            
            obj.weightMap = obj.ydata~=0;
            if obj.ron
                obj.weightMap = obj.weightMap./sqrt(max(obj.ydata,0) + obj.ron^2);           
            end
            obj.ydata = obj.ydata.*obj.weightMap/obj.normFactor;
           

            %4\ Initial guess and bounds        
            obj = fittingSetup(obj,'aoinit',inputs.Results.aoinit,'aobounds',inputs.Results.aobounds,...
                'fitR0',inputs.Results.fitR0,'fitCn2',inputs.Results.fitCn2,...
                'fitGains',inputs.Results.fitGains,'jZernGain',inputs.Results.jZernGain,...
                'statModesFunction',inputs.Results.statModesFunction,'fitStatModes',inputs.Results.fitStatModes,...
                'fitBg',inputs.Results.fitBg,'umax',inputs.Results.umax);
            
            %6\ Define options
            obj.fitOption = optimoptions(@lsqcurvefit,'MaxIter',MaxIter,'TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals,...
                'InitDamping',InitDamping,'Display',display,'SpecifyObjectiveGradient',flagJacobian,'Algorithm',algo);
            
            %7\ Model definition                                   
            if obj.psfr.psf(obj.idSrc).Samp >=1
                nX = obj.psfr.otf.nOtf;
            else
                nX = floor(obj.fov_fit/obj.psfr.psf(obj.idSrc).Samp);
            end
            % Define the fft phasor
            Nx           = obj.psfr.otf.nOtf;
            x            = 1:Nx;
            [X,Y]        = meshgrid(x,x);
            X            = -(X - (Nx/2+1))  * pi*1i/Nx;
            Y            = -(Y - (Nx/2+1))  * pi*1i/Nx;            
            obj.phasor   = @(dx,dy) exp((Y*dy+X*dx));
            obj.modelFUN = @(x,xdata) imageModel(x,xdata,obj);
                      
            %8\ Non-linear least-squares minimization          
            if obj.flagError
                return
            end
            
            if strcmpi(algo,'LEVENBERG-MARQUARDT')
                [obj.beta_ ,~,obj.residualMap_,~,~,~,obj.JacobianMatrix_] = lsqcurvefit(obj.modelFUN,obj.x_init,obj.xdata,obj.ydata,[],[],obj.fitOption);
            else
                [obj.beta_ ,~,obj.residualMap_,~,~,~,obj.JacobianMatrix_] = lsqcurvefit(obj.modelFUN,obj.x_init,obj.xdata,obj.ydata,obj.lbounds,obj.ubounds,obj.fitOption);
            end
        
            %9\ Unpacking results + uncertainties
            obj = updateResults(obj);
                        
        end
    end
    
end
