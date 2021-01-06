function hdr  = defineHeaderFromResults(obj)
inputs = inputParser;
inputs.addRequired('obj',@(x) isa(x,'telemetry') | isa(x,'psfReconstruction') | isa(x,'prime'));
inputs.parse(obj);

if isa(obj,'telemetry')
    trs = obj;
elseif isa(obj,'psfReconstruction')
    trs = obj.trs;
    psfr = obj;
elseif isa(obj,'prime')
    trs = obj.psfr.trs;
    psfr = obj.psfr;
    psfp = obj;
else
    return
end

%Save system status
hdr.SYSGAIN = trs.holoop.gain;
hdr.SYSFREQ = trs.holoop.freq;
hdr.SYSLAT = trs.holoop.lat;
hdr.SYSGAINT = trs.ttloop.gain;
hdr.SYSFREQT = trs.ttloop.freq;
hdr.SYSLATT = trs.ttloop.lat;
hdr.SYSAIRM = trs.tel.airmass;
hdr.SYSEL = trs.tel.elevation;
hdr.SYSAZ = trs.tel.azimuth;
hdr.SYSHLGS = trs.lgs.height*strcmp(trs.aoMode,'LGS');
hdr.SYSEXP = trs.wfs.nExp/trs.holoop.freq;

%Save MASS DIMM
hdr.MASSR0 = trs.atm.r0;
hdr.MASSCN2 = trs.atm.Cn2;
hdr.MASSALT = trs.atm.heights;

 %Save telemetry processing results
 hdr.TRSR0             = trs.res.seeing.r0; %r0 at 500nm but line of sight
 hdr.TRSDR0            = trs.res.seeing.dr0;
 hdr.TRSL0             = trs.res.seeing.L0;
 hdr.TRSDL0            = trs.res.seeing.dL0;
 hdr.TRSSEEV             = trs.res.seeing.seeing_VonKarman;
 hdr.TRSDSEEV            = trs.res.seeing.dseeing_VonKarman;
 hdr.TRSSEEK         = trs.res.seeing.seeing_Kolmo;
 hdr.TRSDSEEK       = trs.res.seeing.dseeing_Kolmo;
 hdr.TRSMNO    = trs.res.noise.std_ho;
 hdr.TRSMNOTT     = trs.res.noise.std_tt;
 hdr.TRSTIP        = std(trs.tipTilt.slopes(1,:),[],2)*1e9;
 hdr.TRSTILT        = std(trs.tipTilt.slopes(1,:),[],2)*1e9;
 hdr.TRSHO         = sqrt(sum(std(trs.rec.res(trs.dm.pupilMask,:),[],2).^2)/nnz(trs.dm.pupilMask))*1e9;
 hdr.TRSNPH             = trs.wfs.nph;
 hdr.TRSNPHTT             = trs.tipTilt.nph;
 
 % Save Gaussian/Moffat fitting results on the sky image
 hdr.SKYSR = trs.sky.SR;
 hdr.SKYDSR = trs.sky.dSR;
 hdr.SKYDFW = trs.sky.dFWHM; 
 hdr.SKYFWX = trs.sky.FWHMx;
 hdr.SKYFWY = trs.sky.FWHMy;
 hdr.SKYDFW = trs.sky.dFWHM;     
 
 hdr.SKYGSR = trs.sky.catalogs.gaussian.SR;
 hdr.SKYGDSR = trs.sky.catalogs.gaussian.dSR;
 hdr.SKYGFWX = trs.sky.catalogs.gaussian.fwhm_x;
 hdr.SKYGDFWX = trs.sky.catalogs.gaussian.dfwhm_x;
 hdr.SKYGFWY = trs.sky.catalogs.gaussian.fwhm_y;
 hdr.SKYGDFWY = trs.sky.catalogs.gaussian.dfwhm_y; 
 hdr.SKYGX = trs.sky.catalogs.gaussian.x;
 hdr.SKYGY = trs.sky.catalogs.gaussian.y;
 hdr.SKYGDX = trs.sky.catalogs.gaussian.dx;
 hdr.SKYGDY = trs.sky.catalogs.gaussian.dy;
 hdr.SKYGF = trs.sky.catalogs.gaussian.flux;
 hdr.SKYGDF = trs.sky.catalogs.gaussian.dflux;
 hdr.SKYGFVU = trs.sky.catalogs.gaussian.fvu; 
 
 hdr.SKYMSR = trs.sky.catalogs.moffat.SR;
 hdr.SKYMDSR = trs.sky.catalogs.moffat.dSR;
  hdr.SKYMFWX = trs.sky.catalogs.moffat.fwhm_x;
 hdr.SKYMDFWX = trs.sky.catalogs.moffat.dfwhm_x;
 hdr.SKYMFWY = trs.sky.catalogs.moffat.fwhm_y;
 hdr.SKYMDFWY = trs.sky.catalogs.moffat.dfwhm_y;
 hdr.SKYMX = trs.sky.catalogs.moffat.x;
 hdr.SKYMY = trs.sky.catalogs.moffat.y;
 hdr.SKYMDX = trs.sky.catalogs.moffat.dx;
 hdr.SKYMDY = trs.sky.catalogs.moffat.dy;
 hdr.SKYMF = trs.sky.catalogs.moffat.flux;
 hdr.SKYMDF = trs.sky.catalogs.moffat.dflux;
 hdr.SKYMFVU = trs.sky.catalogs.moffat.fvu;
 
 % Get PSF results
 if exist('psfr','var')
     
     % PSF metrics
     hdr.RECSR         = psfr.psf.SR;
     hdr.RECFWX      = psfr.psf.FWHMx;
     hdr.RECFWY      = psfr.psf.FWHMy;
     hdr.RECDFWX      = psfr.psf.dFWHM;          
          
     % Get the error breakdown
     wfe = errorBreakDown(psfr);     
     hdr.RECMAR      = wfe.sr_mar;
     hdr.RECPAR      = wfe.sr_par;
     hdr.RECPAR      = wfe.wfe_tot;
     hdr.RECNCPA    = wfe.wfe_ncpa;
     hdr.RECFIT     = wfe.wfe_fit;
     hdr.RECLAG      = wfe.wfe_lag;
     hdr.RECNOI   = wfe.wfe_noise;
     hdr.RECALIAS   = wfe.wfe_alias;
     hdr.RECANISO   = wfe.wfe_aniso;
     hdr.RECTT  = wfe.wfe_tt;
     hdr.RECNOITT = wfe.wfe_noiseTT;
     
     % PSF fitting
     hdr.RECX = psfr.psf.catalogs.ref.x;
     hdr.RECY = psfr.psf.catalogs.ref.y;
     hdr.RECDX = psfr.psf.catalogs.ref.dx;
     hdr.RECDY = psfr.psf.catalogs.ref.dy;
     hdr.RECF = psfr.psf.catalogs.ref.flux;
     hdr.RECDF = psfr.psf.catalogs.ref.dflux;
     hdr.RECFVU = psfr.psf.catalogs.ref.fvu;
 end
 

  if exist('psfp','var')     
     % PSF metrics     
     hdr.PRISR         = psfp.psf.SR;
     hdr.PRIFWX      = psfp.psf.FWHMx;
     hdr.PRIFWY      = psfp.psf.FWHMy;
     hdr.PRIDFWX      = psfp.psf.dFWHM;    
     
     % Get the error breakdown
     wfe = errorBreakDown(psfp);     
     hdr.PRIMAR      = wfe.sr_mar;
     hdr.PRIPAR      = wfe.sr_par;
     hdr.PRIWFE     = wfe.wfe_tot;
     hdr.PRINCPA    = wfe.wfe_ncpa;
     hdr.PRIFIT     = wfe.wfe_fit;
     hdr.PRILAG      = wfe.wfe_lag;
     hdr.PRINOI   = wfe.wfe_noise;
     hdr.PRIALIAS   = wfe.wfe_alias;
     hdr.PRIANISO   = wfe.wfe_aniso;
     hdr.PRITT  = wfe.wfe_tt;
     hdr.PRINOITT = wfe.wfe_noiseTT;
     
     % PSF fitting
     hdr.PRIX = psfp.catalog_fit.x;
     hdr.PRIY = psfp.catalog_fit.y;
     hdr.PRIDX = psfp.catalog_fit.dx;
     hdr.PRIDY = psfp.catalog_fit.dy;
     hdr.PRIF = psfp.catalog_fit.flux;
     hdr.PRIDF = psfp.catalog_fit.dflux;
     hdr.PRIFVU = psfp.catalog_fit.fvu;
     
     % Retrieved parameters
     hdr = setfield(hdr,'PRIGAO',psfp.gains_fit.gAO);
     hdr = setfield(hdr,'PRIDGAO',psfp.gains_fit.dgAO);       
     hdr.PRIGAOR = psfp.gains_fit.gAO(end);     
     hdr.PRIDGAOR = psfp.gains_fit.dgAO(end);     
     hdr.PRIGTT = psfp.gains_fit.gTT;
     hdr.PRIDGTT = psfp.gains_fit.dgTT;
     hdr.PRIGAL = psfp.gains_fit.gAl;
     hdr.PRIDGAL = psfp.gains_fit.dgAl;
     hdr.PRIR0 = psfp.atm_fit.r0*(0.5e-6/psfp.atm_fit.wavelength)^1.2;
     hdr.PRIDR0 = psfp.atm_fit.dr0*(0.5e-6/psfp.atm_fit.wavelength)^1.2;
     hdr = setfield(hdr,'PRICN2',psfp.atm_fit.Cn2);
     hdr = setfield(hdr,'PRIDCN2',psfp.atm_fit.dCn2);    
     if ~isempty(psfp.map_fit)
         hdr.PRICOEF = psfp.map_fit.coefs;
         hdr.PRIDCOEF = psfp.map_fit.dcoefs;
         hdr.SYSFOC = psfp.map_fit.wfsFocusCalib;
         hdr.PRIFOC = psfp.map_fit.wfsFocusMeas;
         hdr.PRIDFOC = psfp.map_fit.dwfsFocusMeas;
     end
    

 end
  
end
