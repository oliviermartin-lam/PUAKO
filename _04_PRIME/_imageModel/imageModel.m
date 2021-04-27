function [im, J_im] = imageModel(x,xdata,obj)  

%% 1\ Split parameters

xS = x(1:3*obj.psfr.trs.src(obj.idSrc).nObj);
if obj.fitBg
    xPSF = x(1+3*obj.psfr.trs.src(obj.idSrc).nObj:end-1);
    xBg = x(end);
else
    xPSF = x(1+3*obj.psfr.trs.src(obj.idSrc).nObj:end);
end
                          
%% 2\ Grab inputs
%2.1. PSF parameters init
r053Init    = (obj.psfr.trs.res.seeing.r0*(obj.psfr.trs.cam.wavelength/0.5e-6)^1.2)^(-5/3);
fov_fit     = obj.fov_fit;
fov_sky     = obj.psfr.trs.cam.resolution;
samp        = obj.psfr.trs.cam.samp;
sf          = obj.psfr.sf;
otf         = obj.psfr.otf;

%2.2. Fixed parameters
r053    = obj.x_fixed.r053;
Cn2     = obj.x_fixed.Cn2;
gAO     = obj.x_fixed.gAO;
gTT     = obj.x_fixed.gTT;
gAl     = obj.x_fixed.gAl;


%2.3. Fitted parameters indexes
obj.idxCn2      = xdata{1};
obj.idxR0       = xdata{2};
obj.idxDao      = xdata{3};
obj.idxDtt      = xdata{4};
obj.idxDal      = xdata{5};
obj.idxStatModes= xdata{6};

%2.4. Update parameters
if ~isempty(obj.idxR0)
    r053= sum(xPSF(obj.idxR0));
end
if ~isempty(obj.idxCn2)
    Cn2 = xPSF(obj.idxCn2);
    r053 = sum(Cn2);
end
if ~isempty(obj.idxDao)
    gAO = xPSF(obj.idxDao);
end
if ~isempty(obj.idxDal)
    gAl = xPSF(obj.idxDal);
end
if ~isempty(obj.idxDtt)
    gTT = xPSF(obj.idxDtt);
end
if ~isempty(obj.idxStatModes)
    statCoefs = xPSF(obj.idxStatModes);
end

%2.5. Stellar parameters
nSrc = length(xS)/3;
fSrc = xS(1:nSrc);
xSrc = xS(1+nSrc:2*nSrc);
ySrc = xS(1+2*nSrc:3*nSrc);
dX   = [xSrc;ySrc];

%if numel(samp)~=numel(xSrc)
%    samp = samp(1)*ones(1,nSrc);
%end

%% 3\ Get the image model accordingly the input parameters

im  = zeros(fov_fit,fov_fit);
if ~isempty(obj.idxStatModes)
    obj.psfr= computeStaticOpticalTransferFunction(obj.psfr,'statModes',{obj.statModes,statCoefs});
end

%3.1 Phase structure function on-axis
Don     = r053/r053Init*(sf.Dfit+ sf.Dal*gAl) + sum(bsxfun(@times,sf.Dao_z,reshape(gAO,1,1,[])), 3) + sf.Dtt*gTT;%*meter2rad^2;
otfOn   = obj.psfr.otf.otfStat.*exp(-0.5*Don).*otf.otfCCD;
sf.Dani = 0;

for iSrc = 1:obj.psfr.trs.src(obj.idSrc).nObj
    %3.2 Anisoplanatism
    if numel(Cn2) > 1 && obj.psfr.flags.isAniso        
        Cn2     = Cn2/sum(Cn2)*r053; %works regardless we fit Cn2 or not
        sf.Dani = sum(bsxfun(@times, squeeze(sf.Dani_l(:,:,:,iSrc)), reshape(Cn2,1,1,[])), 3);  
        %3.3 OTF Multiplication
        otf.Kani = exp(-0.5*sf.Dani);
        otf.otfShannon = otfOn.*otf.Kani;
    else
        otf.otfShannon = otfOn;
    end
       
    %3.4 PSF calculation
    psf_ij  = puakoTools.otfShannon2psf(otf.otfShannon.*obj.phasor(dX(2,iSrc),dX(1,iSrc)),samp,fov_fit);
    S       = sum(sum(puakoTools.crop(psf_ij,fov_sky)));
    psf_ij  = psf_ij/S;
    
    %3.5 Combine sources
    im  = im + abs(fSrc(iSrc))*psf_ij;
end
               
%3.6. Cropping
if fov_sky ~= fov_fit
    im = puakoTools.crop(im,fov_sky);
end

%3.7 Include background           
if obj.fitBg
    im = im + xBg;
end

% 3.8. Weighting
im = im.*obj.weightMap;

%% 4.\ JACOBIAN CALCULATION: TO BE REVIEW AS IT DOES NOT COMPLY WITH THE EMPIRICAL JACOBIAN !!!!!

if nargout >1
    
    %1\ Jacobian with respect to photometry and astrometry - status:ok
    nX = numel(xS) + numel(xPSF);
    if obj.fitBg
        nX = nX + numel(xBg);
    end
    
    J_im = zeros(fov_sky^2,nX);
    for iSrc = 1:nSrc
        % Photometry
        J_im(:,1)   = im(:)/fSrc;
        % Astrometry
        [u,v]       = freqspace(size(im),'meshgrid');
        phasor      = exp(-1i*pi*(v*dX(1) + u*dX(2)));
        otf_sky     = puakoTools.psf2otf(im);
        otf_sky     = otf_sky/max(otf_sky(:));
        J_im(:,2)   = reshape(real(fftshift(ifft2(fftshift(-1i*pi*v.*otf_sky)))).*obj.weightMap,fov_sky^2,1);
        J_im(:,3)   = reshape(real(fftshift(ifft2(fftshift(-1i*pi*u.*otf_sky)))).*obj.weightMap,fov_sky^2,1);
    end
    
    %2\ Jacobian with respect to PSF parameters - status: not ok
    %These ones are not fully ok
    if ~isempty(obj.idxDal)
        gAl = xPSF(obj.idxDal);
    else
        gAl = 1;
    end
    if ~isempty(obj.idxR0)
        r053= sum(xPSF(obj.idxR0));
    else
        r053 = obj.psfr.trs.atm.r0^(-5/3);
    end
       
    n=1;
    if ~isempty(obj.idxR0)
        otf_1               = puakoTools.interpolateOtf(padarray(-0.5*(sf.Dfit+ sf.Dal*gAl).*otf.otfShannon,round((samp-1)*[fov_fit,fov_fit]/2),'both'),fov_sky);
        J_im(:,3*nSrc+n)    = reshape(real(fftshift(ifft2(fftshift(otf_1.*phasor)))).*obj.weightMap,fov_sky^2,[]);
        n = n+1;
    end
    
    if ~isempty(obj.idxDao)
        for k=1:obj.nGainsHO
            otf_1           = puakoTools.interpolateOtf(padarray(-0.5*sf.Dao_z(:,:,k).*otf.otfShannon,round((samp-1)*[fov_fit,fov_fit]/2),'both'),fov_sky);
            J_im(:,3*nSrc+n)= reshape(real(fftshift(ifft2(fftshift(otf_1.*phasor)))).*obj.weightMap,fov_sky^2,[]);
            n = n+1;
        end
    end
    if ~isempty(obj.idxDtt)
        otf_1               = puakoTools.interpolateOtf(padarray(-0.5*sf.Dtt.*otf.otfShannon,round((samp-1)*[fov_fit,fov_fit]/2),'both'),fov_sky);
        J_im(:,3*nSrc+n)    = reshape(real(fftshift(ifft2(fftshift(otf_1.*phasor)))).*obj.weightMap,fov_sky^2,[]);
        n = n+1;
    end
    if ~isempty(obj.idxDal)
        otf_1               = puakoTools.interpolateOtf(padarray(-0.5*r053*sf.Dal.*otf.otfShannon,round((samp-1)*[fov_fit,fov_fit]/2),'both'),fov_sky);
        J_im(:,3*nSrc+n)    = reshape(real(fftshift(ifft2(fftshift(otf_1.*phasor)))).*obj.weightMap,fov_sky^2,[]);
        n = n+1;
    end
    
    %Jacobian with respect to static modes  - status: not ok
     if ~isempty(obj.idxStatModes)
         P              =  puakoTools.enlargePupil(obj.psfr.trs.tel.pupil,2);
         phi_stat       = puakoTools.enlargePupil(obj.psfr.res.static.full_map ,2);
         
         for i=1:obj.jZernGain
             E          = -i*1e-9*2*pi/obj.psfr.trs.cam.wavelength*P.*exp(1i*phi_stat).*obj.zerModes(:,i);             
             otfStat2   = fftshift(puakoTools.fftCorrel(E,E));
             otf2       = otfStat2.*exp(-0.5*Don).*otf.otfCCD;
             if numel(Cn2) > 1 && obj.psfr.flags.isAniso
                 otf2 = otf2.*exp(-0.5*sf.Dani);
             end
             otf_1 = puakoTools.interpolateOtf(padarray(otf2,round((samp-1)*[fov_fit,fov_fit]/2),'both'),fov_sky);
             
             %E = P.*exp(1i*phi_stat);
             %otfStat = fftshift(puakoTools.fftCorrel(E,E));
             J_im(:,3*nSrc+n)    = reshape(real(fftshift(ifft2(fftshift(otf_1.*phasor)))).*obj.weightMap,fov_sky^2,[]);
             n = n+1;
         end
    end
    
    % Jacobian with respect to background - status:ok
    if obj.fitBg
        J_im(:,3*nSrc+n:end) = obj.weightMap(:);
    end
end                  
        
