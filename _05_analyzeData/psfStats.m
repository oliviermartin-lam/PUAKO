classdef psfStats < handle
    
    properties (SetObservable=true)       
        % system
        pupil;
        D;
        % Image
        image;
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
    end
    
    
    methods
        
        function obj = psfStats(image,pupil,wavelength,Samp,psInMas,varargin)
            
            % ---------  Checking inputs
            inputs = inputParser;
            inputs.addRequired('image', @isnumeric );
            inputs.addRequired('pupil', @isnumeric );            
            inputs.addRequired('wavelength', @isnumeric );
            inputs.addRequired('Samp', @isnumeric );
            inputs.addRequired('psInMas', @isnumeric );      
            inputs.addParameter('flagMoffat',false,@(x) isnumeric(x) || islogical(x));
            inputs.addParameter('flagGaussian',false,@(x) isnumeric(x)|| islogical(x));
            inputs.addParameter('includeDiffraction',true,@(x) isnumeric(x)|| islogical(x));
            inputs.addParameter('im_ref',[],@isnumeric);
            inputs.parse(image,pupil,wavelength,Samp,psInMas,varargin{:});

           
            %1\ Parsing input
            obj.pupil        = inputs.Results.pupil;
            obj.pixelScale= inputs.Results.psInMas;
            obj.image      = inputs.Results.image;
            obj.Samp       = inputs.Results.Samp;
            obj.wavelength = inputs.Results.wavelength;
            obj.flags.fitMoffat      = inputs.Results.flagMoffat;
            obj.flags.fitGaussian    = inputs.Results.flagGaussian;
            obj.flags.includeDiffraction    = inputs.Results.includeDiffraction;
            obj.psfResolution   = size(obj.image,1);
            obj.im_ref = inputs.Results.im_ref;

            obj.D = constants.radian2mas*obj.wavelength/(2*obj.Samp*obj.pixelScale);    
            obj.psfFieldOfView = obj.psfResolution * obj.pixelScale;

            %2\ Diffraction limit PSF
            obj.psfDL = tools.telescopePsf(obj.pupil,obj.Samp*2);
            obj.psfDL = tools.crop(obj.psfDL,obj.psfResolution);
            obj.psfDL = obj.psfDL/sum(obj.psfDL(:));            
            
            %3\ Computing image statistics            
            %3.1. Coordinates in the focal plane in mas
            obj.focalGrid  = getGridCoordinates(obj.psfResolution,obj.psfResolution,obj.psfFieldOfView/2);
        
             %3.2 Flux, background value and read-out noise (ADU)
            [obj.flux,obj.bg,obj.ron] = tools.getFlux(obj.image);
            
            % 3.3 Strehl ratio       
            [obj.SR,obj.dSR] = tools.getStrehl(obj.image,obj.pupil,obj.Samp);
            
            % 3.4 FWHM
            [obj.FWHMx,obj.FWHMy,obj.dFWHM,obj.Ellipticity] = tools.getFWHM(obj.image,obj.pixelScale,8,'contour');
            
            %3.5 Ensquared Energy      
            obj.EncircledEnergy = tools.getEncircledEnergy(obj.image);           
            
            %4\ PSF fitting
            if obj.flags.fitMoffat
                obj = obj.getMoffat();                
            end
            
            if obj.flags.fitGaussian                
                obj = obj.getGaussian();               
            end
                  
            if ~isempty(obj.im_ref)
                [tmp,obj.im_fit] = tools.findStellarParameters(obj.im_ref,obj.image,[1,0,0]);
                obj.im_res = obj.im_fit - obj.im_ref;
                obj.catalogs.ref.x = tmp(1,:);
                obj.catalogs.ref.y = tmp(2,:);
                obj.catalogs.ref.flux = tmp(3,:);
                obj.catalogs.ref.dx = tmp(4,:);
                obj.catalogs.ref.dy = tmp(5,:);
                obj.catalogs.ref.dflux = tmp(6,:);          
                obj.catalogs.ref.fvu = tools.getFVU(obj.im_ref,obj.im_fit);
                [obj.SR,obj.dSR] = tools.getStrehl(obj.im_fit,obj.pupil,obj.Samp);
            end
        end
                     
        function obj = getMoffat(obj,xinit)
                        
            %1\ Define the Grid
            normFactor = sum(obj.image(:));
            im   = obj.image/normFactor;
            tmp = obj.focalGrid;        
            xdata = [];
            xdata{1} = tmp.x2D;
            xdata{2} = tmp.y2D;
            
            %2\ Initial guess
            if ~exist('xinit','var')
                xinit = [max(im(:)),obj.pixelScale*2,obj.pixelScale*2,1,0,0,0,0];
                %Amplitude - alpha_x - alpha_y - beta - theta - dx - dy - background
            end
            %3\ Parameters hard-bounds
            lb = [0,0,0,0,-pi,-5*obj.pixelScale,-5*obj.pixelScale,-5*std(im(:))];
            ub = [max(im(:))*10,obj.psfFieldOfView/2,obj.psfFieldOfView/2,20,pi,5*obj.pixelScale,5*obj.pixelScale,5*std(im(:))];            
            
            %4\ Fitting options            
            opt = optimoptions(@lsqcurvefit,'MaxIter',3e2,...
                'TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',3e2,...
                'InitDamping',1,'Display','iter');
            
            %5\ Define the model
            if obj.flags.includeDiffraction                
                f  = @(x,xdata) tools.convolve(tools.moffat(x(1:end-1),xdata),obj.psfDL)+x(end);
            else
                f  = @(x,xdata) tools.moffat(x(1:end-1),xdata)+x(end);
            end

            %6\ Do the fitting
            [beta,~,R,~,~,~,J] = lsqcurvefit(f,xinit,xdata,im,lb,ub,opt);                        
            beta(1) = beta(1)*normFactor;            
            beta(end) = beta(end)*normFactor;            
            
            obj.MoffatStd      = diff(nlparci(beta,R,'jacobian',J),1,2);
            obj.MoffatParam    = beta;
            obj.MoffatImage    = f(obj.MoffatParam,xdata);
            obj.MoffatResidual = im - obj.MoffatImage;
            
            %7\ Getting FWHM and Ellipticity                        
            obj.catalogs.moffat.fwhm_x      = 2*beta(2)*sqrt(2^(1./beta(4))-1); % FWHM-x
            obj.catalogs.moffat.dfwhm_x      = hypot(obj.MoffatStd(2)*2*sqrt(2^(1./beta(4))-1), obj.MoffatStd(4)*beta(2)*log(2)/sqrt(2^(1./beta(4))-1)/beta(4)^2 ) ;
            obj.catalogs.moffat.fwhm_y      = 2*beta(3)*sqrt(2^(1./beta(4))-1); % FWHM - y
            obj.catalogs.moffat.dfwhm_y      =hypot(obj.MoffatStd(3)*2*sqrt(2^(1./beta(4))-1), obj.MoffatStd(4)*beta(3)*log(2)/sqrt(2^(1./beta(4))-1)/beta(4)^2 ) ;
            obj.catalogs.moffat.aspectRatio = max([beta(2)./beta(3) beta(3)./beta(2)]); % ellipticity        
            
             %8\ Astrometry and photometry (ADU)
            obj.catalogs.moffat.x = beta(6);
            obj.catalogs.moffat.dx = obj.MoffatStd(6);
            obj.catalogs.moffat.y = beta(7);
            obj.catalogs.moffat.dy = obj.MoffatStd(7);
            obj.catalogs.moffat.flux = beta(1);
            obj.catalogs.moffat.dflux = obj.MoffatStd(1);
            obj.catalogs.moffat.fvu = tools.getFVU(obj.image,obj.MoffatImage);
        end
        
        function obj = getGaussian(obj,xinit)
            
            %1\ Define the Grid
            normFactor = sum(obj.image(:));
            im   = obj.image/normFactor;
            tmp = obj.focalGrid;
            xdata = [];
            xdata{1} = tmp.x2D;
            xdata{2} = tmp.y2D;
            
             %2\ Initial guess
             if ~exist('xinit','var')
                 xinit = [max(im(:)),obj.pixelScale*2,obj.pixelScale*2,0,0,0,0];
                 %Amplitude - alpha_x - alpha_y - theta - dx - dy
             end
              
             %3\ Parameters hard-bounds
            lb = [0,0,0,-pi,-5*obj.pixelScale,-5*obj.pixelScale,-5*std(im(:))];
            ub = [max(im(:))*10,obj.psfFieldOfView/2,obj.psfFieldOfView/2,pi,5*obj.pixelScale,5*obj.pixelScale,5*std(im(:))];
            
            %4\ Fitting options            
            opt = optimoptions(@lsqcurvefit,'MaxIter',3e2,...
                'TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',3e2,...
                'InitDamping',1,'Display','iter');
            
            %5\ Define the model
            if obj.flags.includeDiffraction                               
                f = @(x,xdata) max(im(:))*tools.convolve(tools.gaussian(x,xdata),obj.psfDL)+x(end);
            else
                f  = @(x,xdata) tools.gaussian(x,xdata)+x(end);
            end
            
            %6\ Do the fitting
            [beta,~,R,~,~,~,J]   = lsqcurvefit(f,xinit,xdata,im,lb,ub,opt);            
            beta(1) = beta(1)*normFactor;
            beta(end) = beta(end)*normFactor;
            
            obj.GaussianStd      = diff(nlparci(beta,R,'jacobian',J),1,2);
            obj.GaussianParam    = beta;
            obj.GaussianImage    = f(obj.GaussianParam,xdata);
            obj.GaussianResidual = obj.image - obj.GaussianImage;
            
            %7\ Getting FWHM and Ellipticity
            obj.catalogs.gaussian.fwhm_x      = 2*beta(2)*sqrt(2*log(2)); % FWHM-x
            obj.catalogs.gaussian.dfwhm_x      = 2*obj.GaussianStd(2)*sqrt(2*log(2)); 
            obj.catalogs.gaussian.fwhm_y      = 2*beta(3)*sqrt(2*log(2)); % FWHM - y
            obj.catalogs.gaussian.dfwhm_y      = 2*obj.GaussianStd(3)*sqrt(2*log(2)); 
            obj.catalogs.gaussian.aspectRatio = max([beta(2)./beta(3) beta(3)./beta(2)]); % ellipticity        
                        
            %8\ Astrometry and photometry (ADU)
            obj.catalogs.gaussian.x = beta(5);
            obj.catalogs.gaussian.dx = obj.GaussianStd(5);
            obj.catalogs.gaussian.y = beta(6);
            obj.catalogs.gaussian.dy = obj.GaussianStd(6);
            obj.catalogs.gaussian.flux = beta(1);
            obj.catalogs.gaussian.dflux = obj.GaussianStd(1);
            obj.catalogs.gaussian.fvu = tools.getFVU(obj.image,obj.GaussianImage);
        end
                
    end
    
end
