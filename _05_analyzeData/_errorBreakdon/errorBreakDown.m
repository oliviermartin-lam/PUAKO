%{
------------HEADER-----------------
Objective          ::  Compute bi-dimensional phase structure function of the
aliasing phase

INPUT VARS
 r0,L0          :: Respectively the Fried's parameter and outer in meter
  d                :: The sub-aperture size
nActu          :: the \# DM actuators
nT               :: number of multiple of the Dm cut-off frequency for
defining the maximal reconstructed frequency
 
OUTPUT VARS
 sfFit_2D             :: bi-dimensional phase structure function map
(Toeplitz)
 psd               :: bi-dimensional fitting PSD
sfFit_4D              : bi-dimensional phase structure matrix
Created by         :: O. Beltramo-Martin - ONERA/LAM
Creation date      :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function wfe = errorBreakDown(obj,varargin)
inputs = inputParser;
inputs.addRequired('obj',@(x) isa(x,'psfReconstruction') || isa(x,'telemetry') || isa(x,'prime'));
inputs.addParameter('display', false,@islogical);
inputs.parse(obj);
display= inputs.Results.display;

%%
if isa(obj,'prime')
    
    % Get inputs
    wvl = obj.psfr.trs.cam.wavelength;
    otf   = obj.psfr.otf;
    sf     = obj.psfr.sf;
     notf  = size(otf.otfDL,1);
    u       = linspace(-1,1-~mod(notf,2)/notf,notf);
    S       = trapz(u,trapz(u,otf.otfDL));
    r0src = obj.r0_fit;
    r053 = r0src^(-5/3);
    r053Init = obj.psfr.trs.res.seeing.r0^(-5/3);
    D = obj.psfr.trs.tel.Dcircle;
    
    %1\ Get the static error from the static map (NCPA)
    srStatic = trapz(u,trapz(u,otf.otfStat))/S;
    wfeStatic = tools.sr2wfe(srStatic,wvl);
    
    %2\ Get the fitting error by integrating the Fitting PSD
    srFit = trapz(u,trapz(u,otf.otfDL.*exp(-0.5*sf.Dfit*r053/r053Init)))/S;
    wfeFit = tools.sr2wfe(srFit,wvl);
    
    %3\ Get the aliasing error from the covariance matrix
    if ~isempty(obj.idxDal)
        gAl  = obj.xao_fit(obj.idxDal - length(obj.idxR0));
    else
        gAl = 1;
    end
    
    srAl = trapz(u,trapz(u,otf.otfDL.*exp(-0.5*sf.Dal*r053/r053Init*gAl)))/S;
    wfeAl = tools.sr2wfe(srAl,wvl);
    
    %4\ Get the noise errors
    %4.1 HO WFS Noise
    if ~isempty(obj.idxR0) && ~isempty(obj.idxDho)
        gHO = obj.xao_fit(obj.idxDho- length(obj.idxR0));
    else
        gHO = 1;
    end
    
    wfeNoise = 0;
    msk = obj.trs.dm.pupilMask;   
    cte = 1e9*sqrt(obj.psfr.trs.holoop.tf.pn/size(obj.psfr.res.noise.Cn_ho(msk,msk),1));
    cn_tmp = obj.psfr_.res.noise.Cn_ho;
    
    for i=1:obj.nGainsHO-1
        tmp = obj.psfr_.trs.mat.Hz{i}*obj.psfr.res.noise.Cn_ho*obj.psfr.trs.mat.Hz{i}';
        wfeNoise  = wfeNoise + gHO(i)*trace(tmp);
        cn_tmp = cn_tmp - tmp;
    end
    wfeNoise = sqrt(wfeNoise + gHO(end)*trace(cn_tmp(msk,msk)))*cte;
    
    %4.2 TT WFS Noise
    if ~isempty(obj.idxR0) && ~isempty(obj.idxDtt)
        gTT = obj.xao_fit(obj.idxDtt - length(obj.idxR0));
    else
        gTT = 1;
    end
    wfeNoiseTT = sqrt(trace(gTT*obj.psfr.trs.ttloop.tf.pn*obj.psfr_.res.noise.Cn_tt))*1e9;
    
    %5\. AO Bandwidth errors
    srLag    = trapz(u,trapz(u,otf.otfDL.*exp(-0.5*sum(bsxfun(@times,sf.Dho_z,reshape(gHO,1,1,[])), 3))))/S;
    wfeLag = sqrt(tools.sr2wfe(srLag,wvl)^2 - wfeNoise^2) ;
    
    %6.  Residual tip-tilt
    srTT    = trapz(u,trapz(u,otf.otfDL.*exp(-0.5*sf.Dtt*gTT)))/S;
    wfeTT = sqrt(tools.sr2wfe(srTT,wvl)^2 - wfeNoiseTT^2);
    
    %7\  Anisoplanatism
    srAni = real(trapz(u,trapz(u,otf.Kani.*otf.otfDL)))/S;
    wfeAni= tools.sr2wfe(srAni, wvl);
    
    SRcam = obj.psf.SR;
    dSRcam = obj.psf.dSR;
    srMar = srStatic.*srFit.*srAl.*srAni.*srTT.*srLag;
    wfeTot = tools.sr2wfe(srMar,wvl);
    %%
elseif isa(obj,'psfReconstruction')
    
    % Get inputs
    wvl = obj.trs.cam.wavelength;
    otf   = obj.otf;
    sf     = obj.sf;
    notf  = size(otf.otfDL,1);
    u       = linspace(-1,1-~mod(notf,2)/notf,notf);
    S       = trapz(u,trapz(u,otf.otfDL));
    
    r0src = obj.trs.res.seeing.r0*(obj.trs.cam.wavelength/0.5e-6)^1.2;
    D = obj.trs.tel.Dcircle;
    
    %1\ Get the static error from the static map (NCPA)
    srStatic = trapz(u,trapz(u,real(otf.otfStat)))/S;
    wfeStatic = tools.sr2wfe(srStatic,wvl);
    
    %2\ Get the fitting error by integrating the Fitting PSD
    srFit =  trapz(u,trapz(u,otf.otfDL.*exp(-0.5*sf.Dfit)))/S;
    wfeFit = tools.sr2wfe(srFit,wvl);
    
    %3\ Get the aliasing error from the covariance matrix        
    srAl =  trapz(u,trapz(u,otf.otfDL.*exp(-0.5*sf.Dal)))/S;
    wfeAl = tools.sr2wfe(srAl,wvl);
    
    %4\ Get the noise errors
    %4.1 HO WFS Noise   
    msk = obj.trs.dm.pupilMask;   
    wfeNoise = 1e9*sqrt(obj.trs.holoop.tf.pn*trace(obj.trs.res.noise.Cn_ho(msk,msk)/nnz(msk)));
    
    %4.2 TT WFS Noise
    wfeNoiseTT = sqrt(trace(obj.trs.ttloop.tf.pn*obj.trs.res.noise.Cn_tt))*1e9;
    
    %5\. AO Bandwidth errors
    srLag    = trapz(u,trapz(u,otf.otfDL.*exp(-0.5*sf.Dho)))/S;
    wfeLag = sqrt(tools.sr2wfe(srLag,wvl)^2 - wfeNoise^2) ;
    
    %6.  Residual tip-tilt
    srTT    =  trapz(u,trapz(u,otf.otfDL.*exp(-0.5*sf.Dtt)))/S;
    wfeTT = sqrt(tools.sr2wfe(srTT,wvl)^2 - wfeNoiseTT^2);
    
    %7\  Anisoplanatism
    srAni =  trapz(u,trapz(u,otf.Kani.*otf.otfDL))/S;
    wfeAni= tools.sr2wfe(srAni, wvl);
    
    SRcam = obj.sky.SR;
    dSRcam = obj.sky.dSR;   
    srMar = srStatic.*srFit.*srAl.*srAni.*srTT.*srLag;
    wfeTot = tools.sr2wfe(srMar,wvl);
    %%
    
elseif isa(obj,'telemetry')
    wvl = obj.cam.wavelength;

    %1\ Get the static error from the static map (NCPA)
    wfeStatic = std(obj.tel.static_map(obj.tel.static_map~=0));
    srStatic = tools.wfe2sr(wfeStatic,wvl);

    %2\ Get the fitting error 
    r0src = obj.res.seeing.r0*(obj.cam.wavelength/0.5e-6)^1.2;
    D = obj.tel.Dcircle;
    wfeFit = sqrt(0.3*(obj.dm.pitch/r0src)^(5/3))*wvl*1e9/2/pi;
    srFit = tools.wfe2sr(wfeFit,wvl);

    %3\ Get the aliasing error 
    wfeAl = wfeFit/sqrt(3);
    srAl= tools.wfe2sr(wfeAl,wvl);

    %4\ Get the noise errors
    %4.1 HO WFS Noise
    wfeNoise = 1e9*sqrt(obj.holoop.tf.pn*obj.res.noise.varn_ho);
    
    %4.2 TT WFS Noise
    wfeNoiseTT = 1e9*sqrt(obj.ttloop.tf.pn*obj.res.noise.varn_tt);
    
    %5\. AO Bandwidth errors
    msk = obj.dm.pupilMask;
    wfeLag = 1e9*sqrt(sum(std(obj.rec.res(msk(:),:),[],2).^2/nnz(msk)));
    wfeLag = sqrt(wfeLag^2 - wfeNoise^2) ;
    srLag = tools.wfe2sr(wfeLag,wvl);

    %6.  Residual tip-tilt
    wfeTT = 1e9*sqrt(sum(std(obj.tipTilt.slopes,[],2).^2));
    wfeTT = sqrt(wfeTT^2 - wfeNoiseTT^2);
    srTT = tools.wfe2sr(wfeTT,wvl);
    
    %7\  Anisoplanatism (need to implement theoretical equations)
    wfeAni = 0;
    srAni = tools.wfe2sr(wfeAni,wvl);
    
    [SRcam,dSRcam] = tools.getStrehl(obj.cam.frame,obj.tel.pupil,obj.cam.samp);
    wfeTot = sqrt(sum(wfeStatic^2 + wfeFit^2 + wfeLag^2 ...
    + wfeNoise^2 + wfeAni^2 + wfeTT^2 +  wfeNoiseTT^2 + wfeAl^2));
    srMar = tools.wfe2sr(wfeTot,wvl);
else
    fprintf('Sorry, I do not recognize the format of the object you provide as entry\n');
end
    
%8\  Total residual
SRho  = srStatic*srFit.*srAni*srLag*srAl;
srPar = SRho/(1+tools.sr2wfe(srTT,wvl)^2*(2*pi*1e-9/wvl)^2) + (1-SRho)/(1+(D/r0src)^2);

% -----------------  DISPLAY
if display
    fprintf('-------------------------------\n');
    fprintf('Wavelength\t\t%.3g micron\n', wvl*1e6)
    fprintf('Strehl Image\t\t%.3g%s +/- %.3g\t\n',1e2*SRcam,'%',1e2*dSRcam);
    fprintf('Strehl Marechal \t%.3g%s\t\n',1e2*srMar,'%');
    fprintf('Strehl Parenti\t\t%.3g%s\t\n',1e2*srPar,'%');
    fprintf('Wavefront error (nm)\t%.4g\t\n',wfeTot);
    fprintf('-------------------------------\n');
    fprintf('NCPA calibration\t%.4g\n',wfeStatic);
    fprintf('Atmospheric Fitting\t%.4g\n',wfeFit);
    fprintf('-------------------------------\n');
    fprintf('Servo-lag\t\t%.4g\n',wfeLag);
    fprintf('WFS noise\t\t%.4g\n',wfeNoise);
    fprintf('WFS aliasing\t\t%.4g\n',wfeAl);
    fprintf('Anisoplanatism\t%.4g\n',max(wfeAni));
    %fprintf('Focal anisoplanatism.\t%.4g\n',wfeAnisoLGS);
    fprintf('-------------------------------\n');
    fprintf('Tip-tilt bandwidth\t%.4g\n',wfeTT);
    fprintf('Tip-tilt noise\t\t%.4g\n',wfeNoiseTT);
    %fprintf('Tip-tilt anisoplanatism\t%g\n',wfeAnisoTT);
    fprintf('-------------------------------\n');
end

wfe.sr_tot      = srMar;
wfe.wfe_tot     = wfeTot;
wfe.wfe_ncpa    = wfeStatic;
wfe.wfe_fit     = wfeFit;
wfe.wfe_lag      = wfeLag;
wfe.wfe_noise   = wfeNoise;
wfe.wfe_alias   = wfeAl;
wfe.wfe_aniso   = wfeAni;
wfe.wfe_tt  = wfeTT;
wfe.wfe_noiseTT = wfeNoiseTT;
