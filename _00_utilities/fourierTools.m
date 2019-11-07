classdef fourierTools < handle
    
   % Tools facility for Fourier calculation
    
    methods (Static)
        
        function out = fftCorrel(x,y)
            nPts = length(x(:));
            out  = ifft2(fft2(x).*conj(fft2(y)))./nPts;
        end
        
        function out = rolledFFT(x)
            nPts = length(x);
            if mod(nPts,2)==0
                out  = fftshift(fft2(fftshift(x)))./nPts^2;
            else
                out  = ifftshift(fft2(ifftshift(x)))./nPts^2;
            end
        end
        
        function out = convolve(object,PSF)
            out = fftshift(ifft2(fft2(object).*fft2(PSF)));
            out = sum(object(:)).*out./sum(out(:));
        end
        
        function out = telescopeOtf(pupil,overSampling)
            extendedPup  = tools.enlargeOtf(pupil,overSampling);
            out = fftshift(tools.fftCorrel(extendedPup,extendedPup));
        end
        
        function out = telescopePsf(pupil,overSampling,peakLocation)
            if nargin <3
                peakLocation = 'middlePixel';
            end
            nSize = size(pupil,1);
            
            if overSampling >1
                otf = tools.telescopeOtf(pupil,overSampling);
                otf = tools.interpolateOtf(otf,nSize);
                out = tools.otf2psf(otf,peakLocation);
            else
                
                otf = tools.telescopeOtf(pupil,2);
                otf = tools.interpolateOtf(otf,nSize/overSampling);
                out = tools.otf2psf(otf);
                out = tools.interpolateOtf(out,nSize);
            end
            
        end
        
        function out = pupil2psf(pupil,phase,overSampling)
            if nargin <2
                phase = 0.*pupil;
            end
            if nargin < 3
                overSampling = 1;
            end
            otf = tools.pupil2otf(pupil,phase,overSampling);
            out = tools.otf2psf(otf);
        end
        
        function out = pupil2otf(pupil,phase,Samp)
            if nargin <2
                phase = 0.*pupil;
            end
            if nargin < 3
                Samp = 1;
            end
            
            if Samp >=1
                P   = tools.enlargeOtf(pupil,2*Samp);
                phi = tools.enlargeOtf(phase,2*Samp);
                E   =  P.*exp(1i.*phi);
                out = fftshift(tools.fftCorrel(E,E));
            else
                P   = tools.enlargeOtf(pupil,2);
                phi = tools.enlargeOtf(phase,2);
                E   =  P.*exp(1i.*phi);
                otf = fftshift(tools.fftCorrel(E,E));
                otf = otf/max(otf(:));
                psf  = real(tools.otfShannon2psf(otf,Samp,size(otf,1)));
                out = tools.psf2otf(psf);
            end
            out = out/max(out(:));
        end
        
        function out = psd2cov(psd,pixelScale)            
            if mod(size(psd,1),2)
                out = fftshift(fft2( fftshift(psd)))*pixelScale^2;
            else
                out = fftshift(fft2( ifftshift(psd)))*pixelScale^2;
            end
        end
        
        function out  = cov2sf(cov)
            out = 2*max(cov(:)) - cov - conj(cov);
        end
        
        function out = sf2otf(sf)
            out  = exp(-0.5 * sf);
        end
        
        function out = otf2psf(otf,peakLocation)
            if nargin <2
                peakLocation = 'middlePixel';
            end
            
            nSize     = size(otf,1);
            [u,v]     = freqspace(nSize,'meshgrid');
            
            if strcmp(peakLocation,'middlePixel')
                % ensures to put the psf max on the middle pixel
                fftPhasor = 1;
            else
                fftPhasor = exp(1i.*pi.*(u+v).*0.5);
                %Max spead out on the four middle pixels of the image
            end
            if  mod(nSize,2) == 0
                out       = real(fftshift(ifft2(fftshift(otf.*fftPhasor))));
            else
                out       = real(fftshift(ifft2(ifftshift(otf.*fftPhasor))));
            end
            out       = out./sum(out(:));
        end
        
        function out = otfShannon2psf(otfShannon,Samp,fov)
            inputs = inputParser;
            inputs.addRequired('otfShannon', @isnumeric);
            inputs.addRequired('Samp', @isnumeric);
            inputs.addRequired('fov', @isnumeric);
            inputs.parse(otfShannon,Samp,fov);
            
            if Samp > 1
                %     % Zero-pad the OTF to get the good PSF pixel scale
                otf    = padarray(otfShannon,round((Samp-1)*size(otfShannon)/2),'both');
                % Interpolate the OTF to crop the PSF to the desired FOV
                otf    = tools.interpolateOtf(otf,fov);
                otf    = otf/max(otf(:));
                out    = tools.otf2psf(otf);
            elseif Samp ==1
                otf    = tools.interpolateOtf(otfShannon,fov);
                otf    = otf/max(otf(:));
                out    = tools.otf2psf(otf);
            else
                % OTF are derived for a Nyquist-sampled PSF
                otf        = tools.interpolateOtf(otfShannon,fov/Samp);
                otf(otf<0) = 0;
                otf        = otf/max(otf(:));
                out        = tools.otf2psf(otf);
                out        = tools.interpolateOtf(out,fov);
            end
        end
        
        function out = psf2otf(psf)
            out = (tools.rolledFFT(psf))/sum(psf(:));
        end
        
        function out = psd2otf(psd,pixelScale)
            cov = tools.psd2cov(psd,pixelScale);
            sf  = tools.cov2sf(cov);
            out = tools.sf2otf(sf);
        end
        
        function out = psd2psf(psd,pixelScale)
            otf = fftshift(tools.psd2otf(psd,pixelScale));
            out = tools.otf2psf(otf);
        end
        
    end
end
