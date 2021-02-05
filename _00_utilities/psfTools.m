classdef psfTools < handle
    % Tools facility for psf characterization
    
    methods (Static)
        
        function out = sr2wfe(SR,lambda)
            out = 1e9*sqrt(-log(SR))*lambda/2/pi;
        end
        
        function out = wfe2sr(wfe,lambda)
            out = exp(-(2*pi*wfe*1e-9/lambda)^2);
        end
        
        function [Flux,bg,ron,msk] = getFlux(psf)
            %Define the inner ellipse
            npx      = size(psf);
            y        = linspace(-1,1,npx(1));
            x        = linspace(-1,1,npx(2));
            [X,Y]    = meshgrid(x,y);
            r        = hypot(X,Y);
            msk      = r>1;
            % Computing the residual background
            psfNoise = psf .*msk;
            bg       = median(psfNoise(msk));
            % Computing the read-out noise
            ron      = std(psfNoise(msk));
            %Computing the normalized flux            
            Flux     = sum(psf(:) - bg);
        end
        
        function [FWHMx,FWHMy,dFWHM,imfit,beta] = getFWHM(psf,varargin)
            
            inputs = inputParser;
            inputs.addRequired('psf',@isnumeric);
            inputs.addParameter('pixelScale',1,@isnumeric);
            inputs.addParameter('rebin',4,@isnumeric);
            inputs.addParameter('method','contour',@ischar);
            inputs.addParameter('pupil',[],@isnumeric);
            inputs.addParameter('samp',2,@isnumeric);
            inputs.addParameter('nStars',1,@isnumeric);
            inputs.addParameter('initGuess',[],@isnumeric);
            inputs.addParameter('flagSymetric',false,@islogical);
            inputs.addParameter('iMax',Inf,@isnumeric);
            inputs.parse(psf,varargin{:});

            pixelScale  = inputs.Results.pixelScale;
            rebin       = inputs.Results.rebin;
            method      = inputs.Results.method;
            pupil       = inputs.Results.pupil;
            samp        = inputs.Results.samp;
            nStars      = inputs.Results.nStars;
            initGuess   = inputs.Results.initGuess;
            flagSymetric= inputs.Results.flagSymetric;
            iMax        = inputs.Results.iMax;
            
            % Gaussian and Moffat fitting are not really efficient on
            % anisoplanatic PSF. Prefer the coutour function in such a
            % case. The cutting method is not compliant to PSF not oriented
            % along x or y-axis.
          
            %Interpolation   
            if rebin ~=1
                im2 = puakoTools.interpolateOtf(psf,rebin*size(psf,1));
            else
                im2 = psf;
            end
            
            if strcmp(method,'cutting')
                % Brutal approach when the PSF is centered and aligned
                % x-axis FWHM
                imx     = im2(:,floor(end/2+1));
                idx     = imx >= (max(imx(:))/2.);
                w       = find(idx==1);
                FWHMx   = (max(w) - min(w))/rebin*pixelScale;%sqrt(4.*sum(idx)/pi)*pixelScale/rebin;
                % y-axis FWHM
                imy     = im2(floor(end/2+1),:);
                idx     = imy >= (max(imy(:))/2.);
                w       = find(idx==1);
                FWHMy   = (max(w) - min(w))/rebin*pixelScale;%sqrt(4.*sum(idx)/pi)*pixelScale/rebin;
                theta   = 0;
                imfit = [];
            elseif strcmp(method,'contour')
                % Contour approach~: something wrong about the ellipse
                % orientation
                C       = contourc(im2,max(im2(:))*[0.5 0.5]);
                if ~isempty(C)
                    % centering the ellispe
                    mx      = [max(C(1,2:end)),max(C(2,2:end))];
                    mn      = [min(C(1,2:end)),min(C(2,2:end))];
                    cent    = (mx+mn)/2;
                    wx      = C(1,2:end) - cent(1);
                    wy      = C(2,2:end) - cent(2);
                    % Get the module
                    r       = hypot(wx,wy)/rebin*pixelScale;
                    % Getting the FWHM
                    FWHMx   = 2*max(r);
                    FWHMy   = 2*min(r);
                    % Getting the ellipse orientation
                    xm      = wx(r == max(r));
                    ym      = wy(r == max(r));
                    theta   = mean(abs(cart2pol(xm,ym)*180/pi));%mean(180*atan(ym./xm)/pi);
                    % Angle are counted positively in the reverse clockwise
                    imfit = [];
                else
                    FWHMx = 0;
                    FWHMy = 0;
                    imfit = [];
                end
                % direction.
            elseif strcmp(method,'Gaussian')
                
                x       = 1:rebin*size(psf,1);
                y       = 1:rebin*size(psf,2);
                x       = x - x(floor(end/2+1));
                y       = y - y(floor(end/2+1));
                [X,Y]   = meshgrid(y,x);
                xdata   = {X,Y};
                if isempty(initGuess)
                    initGuess  = [max(psf(:))*ones(1,nStars)/nStars,2,2,0,zeros(1,2*nStars),0];
                end
                
                if isempty(pupil)
                    if nStars > 1
                        f = @(x,xdata) min(puakoTools.multiAnalyticIsoplanatic(x(1:end-1),xdata,'gaussian',flagSymetric) + x(end),iMax);
                    else
                        f = @(x,xdata) min(puakoTools.gaussian(x(1:end-1),xdata,flagSymetric) + x(end),iMax);
                    end
                else
                    psftel = puakoTools.pupil2psf(pupil,0*pupil,samp);
                    if nStars > 1
                        f = @(x,xdata) min(puakoTools.convolve(psftel,puakoTools.multiAnalyticIsoplanatic(x(1:end-1),xdata,'gaussian',flagSymetric)) + x(end),iMax);
                    else
                        f = @(x,xdata) min(puakoTools.convolve(psftel,puakoTools.gaussian(x(1:end-1),xdata,flagSymetric)) + x(end),iMax);
                    end
                end
                beta    = lsqcurvefit(f,initGuess,xdata,min(im2,iMax));
                FWHMx   = beta(1+nStars)*2*sqrt(2*log(2));
                FWHMy   = beta(2+nStars)*2*sqrt(2*log(2));
                theta   = beta(3+nStars);
                if rebin ~=1
                    x       = 1:size(psf,1);
                    x       = x - x(floor(end/2+1));
                    [X,Y]   = meshgrid(x,x);
                    if ~isempty(pupil)                     
                        psftel = puakoTools.crop(puakoTools.pupil2psf(pupil,0*pupil,samp),size(psf,1));
                    end                    
                    imfit   = f([beta(1:nStars),beta(1+nStars:2+nStars)/rebin,beta(1+2*nStars),beta(2+2*nStars:1+3*nStars)/rebin,0],{X,Y});
                else
                    imfit = f(beta,xdata);
                end
                
            elseif strcmp(method,'Moffat')
                x       = 1:rebin*size(psf,1);
                y       = 1:rebin*size(psf,2);
                x       = x - x(floor(end/2+1));
                y       = y - y(floor(end/2+1));
                [X,Y]   = meshgrid(y,x);
                xdata   = {X,Y};
                                
                if flagSymetric
                    if isempty(initGuess)
                        initGuess  = [max(psf(:))*ones(1,nStars)/nStars,2,2,zeros(1,2*nStars),0];
                    end
                    lb         = [zeros(1,nStars),0,0,-rebin*max(size(psf))*ones(1,2*nStars),-3*max(psf(:))];
                    ub         = [1e10*ones(1,nStars),1e2,1e2,rebin*max(size(psf))*ones(1,2*nStars),3*max(psf(:))];
                else
                    if isempty(initGuess)
                        initGuess  = [max(psf(:))*ones(1,nStars)/nStars,2,2,2,0,zeros(1,2*nStars),0];
                    end
                    lb         = [zeros(1,nStars),0,0,0,-pi,-rebin*max(size(psf))*ones(1,2*nStars),-3*max(psf(:))];
                    ub         = [1e10*ones(1,nStars),1e5,1e5,1e5,pi,rebin*max(size(psf))*ones(1,2*nStars),3*max(psf(:))];
                end
                    
                fitOption = optimoptions(@lsqcurvefit,'MaxIter',5e2,'TolX',1e-14,'TolFun',1e-14,'MaxFunEvals',1e3);
                
                if isempty(pupil)
                    if nStars > 1
                        f   = @(x,xdata) min(puakoTools.multiAnalyticIsoplanatic(x(1:end-1),xdata,'moffat',flagSymetric) + x(end),iMax);
                    else
                        f = @(x,xdata) min(puakoTools.moffat(x(1:end-1),xdata,flagSymetric) + x(end),iMax);
                    end
                else
                    psftel = puakoTools.crop(puakoTools.pupil2psf(pupil,0*pupil,samp),rebin*size(psf,1));
                    if nStars > 1
                        f = @(x,xdata) min(puakoTools.convolve(psftel,puakoTools.multiAnalyticIsoplanatic(x(1:end-1),xdata,'moffat',flagSymetric)) + x(end),iMax);
                    else
                        f = @(x,xdata) min(puakoTools.convolve(psftel,puakoTools.moffat(x(1:end-1),xdata,flagSymetric)) + x(end),iMax);
                    end
                end
                
                % fitting
                beta    = lsqcurvefit(f,initGuess,xdata,min(im2,iMax),lb,ub,fitOption);
                if flagSymetric
                    FWHMx   = 2*beta(1+nStars)*sqrt(2^(1./beta(2+nStars))-1);
                    FWHMy   = 2*beta(1+nStars)*sqrt(2^(1./beta(2+nStars))-1);
                    theta   = beta(3+nStars);
                else
                    FWHMx   = 2*beta(1+nStars)*sqrt(2^(1./beta(3+nStars))-1);
                    FWHMy   = 2*beta(2+nStars)*sqrt(2^(1./beta(3+nStars))-1);
                    theta   = beta(4+nStars);
                end
                
                % image interpolation
                if rebin ~=1
                    x       = 1:size(psf,1);
                    x       = x - x(floor(end/2+1));
                    [X,Y]   = meshgrid(x,x);
                    if ~isempty(pupil)                     
                        psftel = puakoTools.crop(puakoTools.pupil2psf(pupil,0*pupil,samp),size(psf,1));
                    end
                    
                    if flagSymetric
                        imfit   = f([beta(1:nStars),beta(1+nStars)/rebin,beta(2+nStars:3+nStars),beta(4+nStars:end-1)/rebin,0],{X,Y});                        
                    else
                        imfit   = f([beta(1:nStars),beta(1+nStars:2+nStars)/rebin,beta(3+nStars:4+nStars),beta(5+nStars:end-1)/rebin,0],{X,Y});
                    end
                else
                    imfit = f(beta,xdata);
                end
            end
            
            % Get Ellipticity
            dFWHM  = sqrt(2)*pixelScale/rebin/2;
            aRatio = max(FWHMx/FWHMy,FWHMy/FWHMx);
        end
        
        function [SR,dSR] = getStrehl(psf0,pupil,Samp)
            
           
            %1\ Get the Image OTF            
            psf     = puakoTools.recenterPSF(psf0,4);
            otf     = puakoTools.psf2otf(psf);
            otf     = otf/max(otf(:));
           
            %2\ Get the Diffraction-limit OTF
            notf    = size(psf0);
            if Samp >= 1
                otfDL   = puakoTools.telescopeOtf(pupil,2*Samp);
                otfDL   = puakoTools.interpolate(otfDL,[notf(2),notf(1)],'spline');
            else
                otfDL   = puakoTools.telescopeOtf(pupil,2);
                psfDL  = puakoTools.otfShannon2psf(otfDL,2*Samp,notf);
                otfDL  = puakoTools.psf2otf(psfDL);
            end
            otfDL   = otfDL/max(otfDL(:));
            
            
            %3\ Get the Strehl
            %-1+~mod(notf,2)/2/notf:2/notf:1-1/2/notf;
            u   = linspace(-1,1-~mod(notf(1),2)/notf(1),notf(1));
            v   = linspace(-1,1-~mod(notf(2),2)/notf(2),notf(2));
            Mdl = trapz(v,trapz(u,otfDL));
            SR  = abs(trapz(v,trapz(u,otf)))/Mdl;
            
            %4\ Get the uncertainty from the precision on the maximum intensity value and the image sum         
            % note that it does not include the uncertainty du to dark subtraction
            [~,~,ron] = puakoTools.getFlux(psf);                        
            Fim = sum(psf(:));
            Mim = max(psf(:));
	    
            dM = sqrt(ron^2 + Mim);
            dF = sqrt(notf(1)*notf(2)*ron^2 + Fim);
            dSR = 5*SR*(dM/Mim + dF/Fim); % precision
        end
        
        function ee = getEncircledEnergy(psf)                                  
            [~,~,pr] = radialCumulative(psf);
            ee  =cumsum(pr)/sum(psf(:));                                        
        end
        
        function out = getFVU(xtrue,xest,nbox)
            if nargin > 2
                n   = length(xtrue);
                idx = floor(n/2+1-nbox/2):floor(n/2+nbox/2);
                xest = xest(idx,idx);
                xtrue= xtrue(idx,idx);
            end
            MSE = sum(sum((xest-xtrue).^2));
            VarX= sum(sum((xtrue - mean(xtrue(:))).^2));
            out = MSE/VarX;
        end
        
    end
end
