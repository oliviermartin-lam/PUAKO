classdef psfStats < handle
    
    properties (SetObservable=true)       
        % system
        pupil;
        D;
        src;
        % Image
        image;
        image_crop;
        im_ref;
        psfResolution;
        pixelScale;
        Samp;
        wavelength;
        psfFieldOfView;               
        focalGrid;
        % Noise
        flux;
        bg;
        ron;
        % PSF Estimates
        SR;
        dSR;
        FWHMx;
        FWHMy;
        dFWHM;        
        Ellipticity;
        EncircledEnergy;            
        catalogs;
        im_fit;
        im_res;
        umax;
        %Moffat
        MoffatImage;
        MoffatResidual;
        MoffatParam;
        MoffatStats;
        MoffatStd;
        %Gaussian
        GaussianImage;
        GaussianResidual;
        GaussianParam;
        GaussianStats;
        GaussianStd;
        psfDL;
        % flags;
        flags = [];
        fitbg;
    end
    
    
    methods
        
        function obj = psfStats(image,pupil,wavelength,Samp,psInMas,src,varargin)
            
            % ---------  Checking inputs
            inputs = inputParser;
            inputs.addRequired('image', @isnumeric );
            inputs.addRequired('pupil', @isnumeric );            
            inputs.addRequired('wavelength', @isnumeric );
            inputs.addRequired('Samp', @isnumeric );
            inputs.addRequired('psInMas', @isnumeric );      
            inputs.addRequired('src', @isstruct );      
            inputs.addParameter('psfResolution',size(image,1),@isnumeric);            
            inputs.addParameter('flagMoffat',false,@(x) isnumeric(x) || islogical(x));
            inputs.addParameter('flagGaussian',false,@(x) isnumeric(x)|| islogical(x));
            inputs.addParameter('includeDiffraction',false,@(x) isnumeric(x)|| islogical(x));
            inputs.addParameter('im_ref',[],@isnumeric);
            inputs.addParameter('fitbg',false,@islogical);
            inputs.addParameter('ron',0,@isnumeric);
            inputs.addParameter('umax',5,@isnumeric);
            inputs.addParameter('x0',[],@isnumeric);
            inputs.parse(image,pupil,wavelength,Samp,psInMas,src,varargin{:});

           
            %1\ Parsing input
            obj.pupil                   = inputs.Results.pupil;
            obj.pixelScale              = inputs.Results.psInMas;
            obj.image                   = inputs.Results.image;
            obj.Samp                    = inputs.Results.Samp;
            obj.src                     = inputs.Results.src;
            obj.psfResolution           = inputs.Results.psfResolution;                                    
            obj.wavelength              = inputs.Results.wavelength;
            obj.flags.fitMoffat         = inputs.Results.flagMoffat;
            obj.flags.fitGaussian       = inputs.Results.flagGaussian;
            obj.flags.includeDiffraction= inputs.Results.includeDiffraction;
            obj.im_ref                  = inputs.Results.im_ref;
            obj.fitbg                   = inputs.Results.fitbg;            
            obj.D                       = constants.radian2mas*obj.wavelength/(2*obj.Samp*obj.pixelScale);    
            obj.psfFieldOfView          = obj.psfResolution * obj.pixelScale;
            obj.ron                     = inputs.Results.ron;
            obj.umax                    = inputs.Results.umax;
            x0                          = inputs.Results.x0;
            % cropping
            if obj.psfResolution <size(obj.image,1)
                obj.image_crop = puakoTools.crop(obj.image,obj.psfResolution);
            else
                obj.image_crop = obj.image;
            end
            
            %2\ Diffraction limit PSF
            obj.psfDL = puakoTools.telescopePsf(obj.pupil,obj.Samp*2);
            obj.psfDL = puakoTools.crop(obj.psfDL,obj.psfResolution);
            obj.psfDL = obj.psfDL/sum(obj.psfDL(:));            
            
            %3\ Computing image statistics            
            %3.1. Coordinates in the focal plane in pixels
            obj.focalGrid  = getGridCoordinates(obj.psfResolution,obj.psfResolution,obj.psfResolution/2);
        
             %3.2 Flux, background value and read-out noise (ADU)
            [obj.flux,obj.bg] = puakoTools.getFlux(obj.image_crop);
            
            % 3.3 Strehl ratio       
            [obj.SR,obj.dSR] = puakoTools.getStrehl(obj.image_crop,obj.pupil,obj.Samp);
            
            % 3.4 FWHM
            [obj.FWHMx,obj.FWHMy,obj.dFWHM,obj.Ellipticity] = puakoTools.getFWHM(obj.image_crop,'pixelScale',obj.pixelScale,'rebin',4,'method','contour');
            
            %3.5 Ensquared Energy      
            obj.EncircledEnergy = puakoTools.getEncircledEnergy(obj.image_crop);           
            
            %4\ PSF fitting
            if obj.flags.fitMoffat
                obj = obj.getMoffat(x0);                
            end
            
            if obj.flags.fitGaussian                
                obj = obj.getGaussian(x0);               
            end
                  
            if ~isempty(obj.im_ref)               
                xS = zeros(1,obj.src.nObj);
                yS = zeros(1,obj.src.nObj);
                fS = ones(1,obj.src.nObj);
                xinit = [xS,yS,fS];
                if obj.fitbg
                    [tmp,obj.im_fit] = puakoTools.findStellarParameters(obj.im_ref,obj.image,[xinit,0]);
                else
                    [tmp,obj.im_fit] = puakoTools.findStellarParameters(obj.im_ref,obj.image,xinit);
                end
                
                if ~isempty(obj.im_fit)
                    obj.im_res = obj.im_fit - obj.im_ref;
                    obj.catalogs.ref.x = tmp(1,:)*obj.pixelScale;
                    obj.catalogs.ref.y = tmp(2,:)*obj.pixelScale;
                    obj.catalogs.ref.flux = tmp(3,:);
                    obj.catalogs.ref.dx = tmp(4,:)*obj.pixelScale;
                    obj.catalogs.ref.dy = tmp(5,:)*obj.pixelScale;
                    obj.catalogs.ref.dflux = tmp(6,:);
                    obj.catalogs.ref.fvu = puakoTools.getFVU(obj.im_ref,obj.im_fit);
                    if obj.fitbg
                        obj.catalogs.ref.bg = tmp(7,:);
                    end
                    if sum(obj.im_fit(:))~=0
                        [obj.SR,obj.dSR] = puakoTools.getStrehl(obj.im_fit,obj.pupil,obj.Samp);
                    end
                end
            end
        end
                     
        function obj = getMoffat(obj,xinit)
                        
            %1\ Define the Grid
            ydata       = obj.image;
            normFactor  = sum(ydata(ydata>0));
            weightMap   = ydata~=0;
            if obj.ron
                weightMap = weightMap./sqrt(max(ydata,0) + obj.ron^2);
            end
            ydata = ydata.*weightMap/normFactor;
            
            %normFactor  = sum(obj.image(:));
            %im          = obj.image/normFactor;
            tmp         = obj.focalGrid;
            xdata       = [];
            xdata{1}    = tmp.y2D;
            xdata{2}    = tmp.x2D;
            nS          = numel(obj.src.x);
      
            %2\ Initial guess
            
            if nargin<1
                xS          = 1e3*(obj.src.x - obj.src.x(end))/obj.pixelScale;
                yS          = 1e3*(obj.src.y - obj.src.y(end))/obj.pixelScale;
                [yref,xref] = find(ydata == max(ydata(:)));
                xS          = xS - (obj.psfResolution/2+1 - xref*ones(1,nS));
                yS          = yS - (obj.psfResolution/2+1 - yref*ones(1,nS));                
                xinit = [max(ydata(:))/nS*ones(1,nS),2,2,1,0,yS,xS,0];
                %Amplitude - alpha_x - alpha_y - beta - theta - dx - dy - background            else
                
                xS = xinit(4+2*nS:4+3*nS-1);
                yS = xinit(4+nS:4+2*nS-1);
            end
            
            %3\ Parameters hard-bounds
            lb = [zeros(1,nS),0,0,0,-pi,yS - obj.umax,xS - obj.umax,-1e10];
            ub = [max(ydata(:))*10*ones(1,nS),obj.psfResolution/2,obj.psfResolution/2,pi,20,yS + obj.umax,xS + obj.umax,1e10];
            
            %4\ Fitting options            
            opt = optimoptions(@lsqcurvefit,'MaxIter',3e2,...
                'TolX',1e-15,'TolFun',1e-15,'MaxFunEvals',1e3,...
                'InitDamping',1,'Display','iter');
            
            %5\ Define the model
            if obj.flags.includeDiffraction                
                f  = @(x,xdata) weightMap.*(puakoTools.convolve(puakoTools.multiAnalyticIsoplanatic(x(1:end-1),xdata,'moffat'),obj.psfDL)+x(end));
            else
                f  = @(x,xdata) weightMap.(puakoTools.multiAnalyticIsoplanatic(x(1:end-1),xdata,'moffat')+x(end));
            end

            %6\ Do the fitting
            [beta,~,R,~,~,~,J] = lsqcurvefit(f,xinit,xdata,im,lb,ub,opt);                        
            beta(1:nS)  = beta(1:nS)*normFactor;            
            beta(end)   = beta(end)*normFactor;            
            
            obj.MoffatStd      = diff(nlparci(beta,R,'jacobian',J),1,2)';
            obj.MoffatParam    = beta;
            weightMap          = 1;
            obj.MoffatImage    = f(obj.MoffatParam,xdata);
            obj.MoffatResidual = im - obj.MoffatImage;
            
            %7\ Getting FWHM and Ellipticity                        
            obj.catalogs.moffat.fwhm_x          = 2*beta(1+nS)*sqrt(2^(1./beta(3+nS))-1)*obj.pixelScale; % FWHM-x
            obj.catalogs.moffat.dfwhm_x         = hypot(obj.MoffatStd(1+nS)*2*sqrt(2^(1./beta(3+nS))-1), obj.MoffatStd(3+nS)*beta(1+nS)*log(2)/sqrt(2^(1./beta(3+nS))-1)/beta(3+nS)^2 ) *obj.pixelScale;
            obj.catalogs.moffat.fwhm_y          = 2*beta(2+nS)*sqrt(2^(1./beta(3+nS))-1)*obj.pixelScale; % FWHM - y
            obj.catalogs.moffat.dfwhm_y         = hypot(obj.MoffatStd(2+nS)*2*sqrt(2^(1./beta(3+nS))-1), obj.MoffatStd(3+nS)*beta(2+nS)*log(2)/sqrt(2^(1./beta(3+nS))-1)/beta(3+nS)^2 ) *obj.pixelScale;
            obj.catalogs.moffat.aspectRatio     = max([beta(1+nS)./beta(2+nS) beta(2+nS)./beta(1+nS)]); % ellipticity        
            
            %8\ Get the Strehl
            [obj.catalogs.moffat.SR,obj.catalogs.moffat.dSR] = puakoTools.getStrehl(obj.MoffatImage,obj.pupil,obj.Samp);
            
             %9\ Astrometry and photometry (ADU)
            obj.catalogs.moffat.x       = beta(5+nS:5+2*nS-1)*obj.pixelScale;
            obj.catalogs.moffat.dx      = obj.MoffatStd(5+nS:5+2*nS-1)*obj.pixelScale;
            obj.catalogs.moffat.y       = beta(5+2*nS:5+3*nS-1)*obj.pixelScale;
            obj.catalogs.moffat.dy      = obj.MoffatStd(5+2*nS:5+3*nS-1)*obj.pixelScale;
            obj.catalogs.moffat.flux    = beta(1:nS);
            obj.catalogs.moffat.dflux   = obj.MoffatStd(1:nS);
            obj.catalogs.moffat.fvu     = puakoTools.getFVU(obj.image,obj.MoffatImage);
        end
        
        function obj = getGaussian(obj,xinit)
            
            %1\ Define the Grid
            ydata       = obj.image;
            normFactor  = sum(ydata(ydata>0));
            weightMap   = ydata~=0;
            if obj.ron
                weightMap = weightMap./sqrt(max(ydata,0) + obj.ron^2);
            end
            ydata = ydata.*weightMap/normFactor;
            
            %normFactor  = sum(obj.image(:));
            %im          = obj.image/normFactor;
            tmp         = obj.focalGrid;
            xdata       = [];
            xdata{1}    = tmp.y2D;
            xdata{2}    = tmp.x2D;
            nS          = numel(obj.src.x);
            % init
           
            %2\ Initial guess
            if nargin<1
                xS          = 1e3*(obj.src.x - obj.src.x(end))/obj.pixelScale;
                yS          = 1e3*(obj.src.y - obj.src.y(end))/obj.pixelScale;
                [yref,xref] = find(ydata == max(ydata(:)));
                xS          = xS - (obj.psfResolution/2+1 - xref*ones(1,nS));
                yS          = yS - (obj.psfResolution/2+1 - yref*ones(1,nS));
                
                xinit = [max(ydata(:))/nS*ones(1,nS),2,2,0,yS,xS,0];
                %Amplitude - alpha_x - alpha_y - theta - dx - dy                
            else
                
                xS = xinit(4+2*nS:4+3*nS-1);
                yS = xinit(4+nS:4+2*nS-1);
            end
            
            %3\ Parameters hard-bounds
            lb = [zeros(1,nS),0,0,-pi,yS - obj.umax,xS - obj.umax,-1e10];
            ub = [1e20*ones(1,nS),obj.psfResolution/2,obj.psfResolution/2,pi,yS + obj.umax,xS + obj.umax,1e10];
            
            %4\ Fitting options
            opt = optimoptions(@lsqcurvefit,'MaxIter',3e2,...
                'TolX',1e-18,'TolFun',1e-18,'MaxFunEvals',1e3,...
                 'Display','iter');
            
            %5\ Define the model
            if obj.flags.includeDiffraction
                f = @(x,xdata) weightMap.*(puakoTools.convolve(puakoTools.multiAnalyticIsoplanatic(x(1:end-1),xdata,'gaussian'),obj.psfDL)+x(end));
            else
                f  = @(x,xdata) weightMap.*(puakoTools.multiAnalyticIsoplanatic(x(1:end-1),xdata,'gaussian')+x(end));
            end
            
            %6\ Do the fitting
            [beta,~,R,~,~,~,J]   = lsqcurvefit(f,xinit,xdata,ydata,lb,ub,opt);
            beta(1:nS) = beta(1:nS)*normFactor;
            beta(end) = beta(end)*normFactor;
            
            obj.GaussianStd      = diff(nlparci(beta,R,'jacobian',J),1,2)';
            obj.GaussianParam    = beta;
            weightMap            = 1;
            obj.GaussianImage    = f(obj.GaussianParam,xdata);
            obj.GaussianResidual = obj.image - obj.GaussianImage;
            
            %7\ Getting FWHM and Ellipticity
            obj.catalogs.gaussian.fwhm_x        = 2*beta(1+nS)*sqrt(2*log(2))*obj.pixelScale; % FWHM-x
            obj.catalogs.gaussian.dfwhm_x       = 2*obj.GaussianStd(1+nS)*sqrt(2*log(2))*obj.pixelScale;
            obj.catalogs.gaussian.fwhm_y        = 2*beta(2+nS)*sqrt(2*log(2))*obj.pixelScale; % FWHM - y
            obj.catalogs.gaussian.dfwhm_y       = 2*obj.GaussianStd(2+nS)*sqrt(2*log(2))*obj.pixelScale;
            obj.catalogs.gaussian.aspectRatio   = max([beta(1+nS)./beta(2+nS) beta(2+nS)./beta(1+nS)]); % ellipticity
            
            %8\ Get the Strehl
            [obj.catalogs.gaussian.SR,obj.catalogs.gaussian.dSR] = puakoTools.getStrehl(obj.GaussianImage,obj.pupil,obj.Samp);
            
            %9\ Astrometry and photometry (ADU)
            obj.catalogs.gaussian.y     = beta(4+nS:4+2*nS-1)*obj.pixelScale;
            obj.catalogs.gaussian.dy    = obj.GaussianStd(4+nS:4+2*nS-1)*obj.pixelScale;
            obj.catalogs.gaussian.x     = beta(4+2*nS:4+3*nS-1)*obj.pixelScale;
            obj.catalogs.gaussian.dx    = obj.GaussianStd(4+2*nS:4+3*nS-1)*obj.pixelScale;
            obj.catalogs.gaussian.flux  = beta(1:nS);
            obj.catalogs.gaussian.dflux = obj.GaussianStd(1:nS);
            obj.catalogs.gaussian.fvu   = puakoTools.getFVU(obj.image,obj.GaussianImage);
        end
                
    end
    
end
