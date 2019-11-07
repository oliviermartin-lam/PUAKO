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
            %Define the inner circle
            npx      = length(psf);
            x        = linspace(-1,1,npx);
            [X,Y]    = meshgrid(x);
            r        = hypot(X,Y);
            msk      = r>1;
            % Computing the residual background
            psfNoise = psf .*msk;
            bg       = median(psfNoise(:));
            % Computing the read-out noise
            ron      = std(psfNoise(:));
            %Computing the normalized flux
            psfFlux  = psf.*(r<=1);
            Flux     = sum(psfFlux(:) -bg);
        end
        
        function [FWHMx,FWHMy,dFWHM,aRatio,theta,beta] = getFWHM(psf,pixelScale,rebin,method)
            
            % Gaussian and Moffat fitting are not really efficient on
            % anisoplanatic PSF. Prefer the coutour function in such a
            % case. The cutting method is not compliant to PSF not oriented
            % along x or y-axis.
            
            if nargin < 3
                rebin = 4;
            end
            if nargin < 4
                method = 'contour';
            end
            %Interpolation   100*tools.getStrehl(camImg,pupMat,overSampling)
            im2     = tools.interpolateOtf(psf,rebin*size(psf,1));
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
                else
                    FWHMx = 0;
                    FWHMy = 0;
                end
                % direction.
            elseif strcmp(method,'Gaussian')
                xdata   = tools.getFocalGrid(size(psf),pixelScale);
                x0      = [max(psf(:)),20,20,1,0,0,0];
                f       = @(x,xdata) tools.gaussianModel(x,xdata);
                beta    = lsqcurvefit(f,x0,xdata,psf);
                FWHMx   = beta(2)*2*sqrt(2*log(2));
                FWHMy   = beta(3)*2*sqrt(2*log(2));
                theta   = beta(4);
            elseif strcmp(method,'Moffat')
                xdata   = tools.getFocalGrid(size(psf),pixelScale);
                x0      = [max(psf(:)),20,20,1,0,0,0];
                f       = @(x,xdata) tools.moffatModel(x,xdata);
                beta    = lsqcurvefit(f,x0,xdata,psf);
                FWHMx   = 2*beta(2)*sqrt(2^(1./beta(4))-1);
                FWHMy   = 2*beta(3)*sqrt(2^(1./beta(4))-1);
                theta   = beta(5);
            end
            % Get Ellipticity
            dFWHM  = sqrt(2)*pixelScale/rebin/2;
            aRatio = max(FWHMx/FWHMy,FWHMy/FWHMx);
        end
        
        function [SR,dSR] = getStrehl(psf0,pupil,Samp)
            
            %1\ Remove the background if any left
            [~,bg] = tools.getFlux(psf0);
            psf0 = psf0 - bg;
            
            %2\ Get the Diffraction-limit OTF
            notf    = size(psf0,2);
            if Samp >= 1
                otfDL   = tools.telescopeOtf(pupil,2*Samp);
                otfDL   = tools.interpolate(otfDL,notf,'spline');
            else
                otfDL   = tools.telescopeOtf(pupil,2);
                psfDL  = tools.otfShannon2psf(otfDL,2*Samp,notf);
                otfDL  = tools.psf2otf(psfDL);
            end
            otfDL   = otfDL/max(otfDL(:));
            
            %3\ Get the Image OTF
            psf     = tools.recenterPSF(psf0,4);
            psf     = psf.*(psf>0);
            F        = sum(psf(:));
            otf     = tools.psf2otf(psf);
            otf     = otf/max(otf(:));

            %4\ Get the Strehl
            %-1+~mod(notf,2)/2/notf:2/notf:1-1/2/notf;
            u       = linspace(-1,1-~mod(notf,2)/notf,notf);
            S       = trapz(u,trapz(u,otfDL));
            M      = abs(trapz(u,trapz(u,otf)));
            SR     = abs(trapz(u,trapz(u,otf)))/S;
            
            %5\ Get the uncertainty from the precision on the maximum intensity value and the image sum         
            % note that it does not include the uncertainty du to dar
            % subtraction
            [~,~,ron] = tools.getFlux(psf);            
            M = M*F;
            dM =sqrt(ron^2 + M)/(S*F);
            dF = M*sqrt(notf^2*ron^2 + F)/F^2/S;
            dSR = 3*hypot(dM,dF); % accuracy 
        end
        
        function ee = getEncircledEnergy(psf)                                  
            [~,~,pr] = radial(psf);
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