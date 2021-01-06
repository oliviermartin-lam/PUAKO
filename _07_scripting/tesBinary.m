clear all;close all;clc;
date = '20170315';

[~,msg] = unix('echo "$USER"');
if contains(msg,'psfr')    
    path_workspace = '/home/psfr/PUAKO/PUAKO/';
    addpath(genpath(path_workspace));
    path_root = '/net/vm-psfr/psfrData/NIRC2/';
    path_calib = '/net/vm-psfr/psfrData/CALIBRATION/';
    path_save='/net/vm-psfr/psfrData/PUAKO_results/';
    path_imag = [path_root,date,'/IMAG/'];   
    path_trs = [path_root,date,'/TRS/'];   
elseif contains(msg,'omartin')
    path_workspace = '/home/omartin/Projects/KVS/CODES/PUAKO/';
    addpath(genpath(path_workspace));
    path_root = '/run/media/omartin/OlivierMartinHDD/DATA/KECK_DATA/';
    path_calib= [path_root,'CALIBRATION/'];
    path_save ='/home/omartin/Projects/KVS/_results/binaryTests/';
    path_imag = [path_root,date,'/IMAG/Binary/'];   
    path_trs  = [path_root,date,'/TRS/'];   
    path_dark = [path_root,date,'/IMAG/Dark/'];
    path_sky  = [path_root,date,'/IMAG/Sky/'];
elseif contains(msg,'sragland')
    path_workspace = '/home/sragland/PUAKO/PUAKO/';
    addpath(genpath(path_workspace));
    path_root = '/net/vm-psfr/psfrData/NIRC2/';
    path_calib = '/net/vm-psfr/psfrData/CALIBRATION/';
    path_save='/net/vm-psfr/psfrData/PUAKO_RESULTS_SAM/RESULTS/';    
    path_imag = [path_root,date,'/IMAG/'];   
    path_trs = [path_root,date,'/TRS/'];   
end
p = puako('path_imag',path_imag,'path_trs',path_trs,'path_calibration',path_calib,'path_dark',path_dark,'path_sky',path_sky);

objname = {p.data_folder_id.imag.name};
nObj = numel(objname);


%% PRIME
iBin        = 1;
SEP_BIN     = zeros(5,nObj);
DMAG_BIN    = zeros(5,nObj);
ANG_BIN     = zeros(5,nObj);
prime_param_piston = zeros(43,nObj);
prime_param_petal = zeros(13,nObj);
nRes        = 100;
im_fit      = zeros(nRes,nRes,6,nObj);
umax        = 5;
ron         = 0;
mode        = 'petal';

for kObj = 1:nObj
   
    t0 = tic();

    % --------------  Running PRIME on the brightest isolated PSF
    
    %1.1. With petal modes
    p.getPrimePSF('objname',objname{kObj},'resolution',150,'fov',170,'ron',300,...
        'fitStatModes',1:6,'statModesFunction','petal','flagBrightestStar',true,...
        'fitGains',[true,true,false],'umax',umax);
    prime_param_petal(:,kObj) = p.psfp.x_final;
    % petal-prime psf
    p.psfp.weightMap = 1;
    psf_prime_petal  = imageModel([1,0,0,p.psfp.x_final(4:end-1),0],p.psfp.xdata,p.psfp);
    psf_prime_petal  = psf_prime_petal/sum(psf_prime_petal(:));
    
    %1.2. With piston modes
    p.getPrimePSF('objname',objname{kObj},'resolution',150,'fov',170,'ron',300,...
        'fitStatModes',1:36,'statModesFunction','piston','flagBrightestStar',true,...
        'fitGains',[true,true,false],'umax',umax);
    prime_param_piston(:,kObj) = p.psfp.x_final;
    % piston-prime psf
    p.psfp.weightMap %initial guess
    idMax                   = find(p.trs.src(iBin).F == max(p.trs.src(iBin).F));
    idMin                   = find(p.trs.src(iBin).F ~= max(p.trs.src(iBin).F));
    pos_init                = -[abs(p.trs.src(iBin).y(idMin) - p.trs.src(iBin).y(idMax)),0,abs(p.trs.src(iBin).x(idMin) - p.trs.src(iBin).x(idMax)),0];
    pos_init                = 1e3*pos_init/p.trs.cam.pixelScale;   
    [yref,xref]             = find(im_bin == max(im_bin(:)));
    pos_init                = (pos_init - (p.trs.cam.resolution/2+1 - [yref*ones(1,2),xref*ones(1,2)]));
    photo_init              = [0.5,0.5];
    bg_init                 = 0;
    psf_prime_piston  = imageModel([1,0,0,p.psfp.x_final(4:end-1),0],p.psfp.xdata,p.psfp);
    psf_prime_piston  = psf_prime_piston/sum(psf_prime_piston(:));
    
    %1.3 Get the PSFR, PRIME-calibrated PSF and Gl 569 A
    psf_psfr   = puakoTools.crop(p.psfr.rec_,p.trs.cam.resolution);
    psf_psfr   = psf_psfr/sum(psf_psfr(:));
    % on-axis extracted PSF
    psf_Gl569A = p.trs.cam.image;
    psf_Gl569A = psf_Gl569A/sum(psf_Gl569A(:));

    % ------------------------- FIT THE BINARY IMAGE 
    %Select the binary image
    p.trs.cam.resolution    = nRes;
    p.trs                   = processDetectorImage(p.trs,'flagGaussian',false,'flagMoffat',false,'ron',ron,'umax',3);
    im_bin                  = p.trs.cam.image(:,:,iBin);
    im_fit(:,:,1,kObj)      = im_bin;
    
    %initial guess
    idMax                   = find(p.trs.src(iBin).F == max(p.trs.src(iBin).F));
    idMin                   = find(p.trs.src(iBin).F ~= max(p.trs.src(iBin).F));
    pos_init                = -[abs(p.trs.src(iBin).y(idMin) - p.trs.src(iBin).y(idMax)),0,abs(p.trs.src(iBin).x(idMin) - p.trs.src(iBin).x(idMax)),0];
    pos_init                = 1e3*pos_init/p.trs.cam.pixelScale;   
    [yref,xref]             = find(im_bin == max(im_bin(:)));
    pos_init                = (pos_init - (p.trs.cam.resolution/2+1 - [yref*ones(1,2),xref*ones(1,2)]));
    photo_init              = [0.5,0.5];
    bg_init                 = 0;

    % Double Gaussian-components fit
    p.trs.sky(iBin).umax    = umax;
    x0                      = [photo_init,3,3,0,pos_init,bg_init];
    p.trs.sky(iBin)         = p.trs.sky(iBin).getGaussian(x0);
    im_fit(:,:,2,kObj)      = p.trs.sky(iBin).GaussianImage;
    SEP_BIN(1,kObj)         = hypot(diff(p.trs.sky(iBin).catalogs.gaussian.x),diff(p.trs.sky(iBin).catalogs.gaussian.y));
    DMAG_BIN(1,kObj)        = 2.5*log10(max(p.trs.sky(iBin).catalogs.gaussian.flux)/min(p.trs.sky(iBin).catalogs.gaussian.flux));
    ANG_BIN(1,kObj)         = 180 - 180/pi*atan(diff(p.trs.sky(iBin).catalogs.gaussian.x)/diff(p.trs.sky(iBin).catalogs.gaussian.y));
  
    % PSF-R   
    [pp,im_fit(:,:,3,kObj),iminit] = puakoTools.findStellarParameters(im_bin,psf_psfr,[pos_init,photo_init,0],'umax',umax,'ron',ron);
    y                       = pp(1,:);
    x                       = pp(2,:);
    flux                    = pp(3,:);
    SEP_BIN(2,kObj)         = hypot(diff(x),diff(y))*p.trs.cam.pixelScale;
    DMAG_BIN(2,kObj)        = 2.5*log10(max(flux)/min(flux));
    ANG_BIN(2,kObj)         = 180 - 180/pi*atan(diff(x)/diff(y));
    
    % PETAL-PRIME PSF
    [pp,im_fit(:,:,4,kObj)] = puakoTools.findStellarParameters(im_bin,psf_prime_petal,[pos_init,photo_init,0],'umax',umax,'ron',ron);
    y                       = pp(1,:);
    x                       = pp(2,:);
    flux                    = pp(3,:);
    SEP_BIN(3,kObj)         = hypot(diff(x),diff(y))*p.trs.cam.pixelScale;
    DMAG_BIN(3,kObj)        = 2.5*log10(max(flux)/min(flux));
    ANG_BIN(3,kObj)         = 180 - 180/pi*atan(diff(x)/diff(y));
    
    % PISTON-PRIME PSF
    [pp,im_fit(:,:,5,kObj)] = puakoTools.findStellarParameters(im_bin,psf_prime_piston,[pos_init,photo_init,0],'umax',umax,'ron',ron);
    y                       = pp(1,:);
    x                       = pp(2,:);
    flux                    = pp(3,:);
    SEP_BIN(4,kObj)         = hypot(diff(x),diff(y))*p.trs.cam.pixelScale;
    DMAG_BIN(4,kObj)        = 2.5*log10(max(flux)/min(flux));
    ANG_BIN(4,kObj)         = 180 - 180/pi*atan(diff(x)/diff(y));
    
    %Gl 569 A          
    [pp,im_fit(:,:,6,kObj)] = puakoTools.findStellarParameters(im_bin,psf_Gl569A,[pos_init,photo_init,0],'umax',umax,'ron',ron);
    y                       = pp(1,:);
    x                       = pp(2,:);
    flux                    = pp(3,:);
    SEP_BIN(5,kObj)         = hypot(diff(x),diff(y))*p.trs.cam.pixelScale;
    DMAG_BIN(5,kObj)        = 2.5*log10(max(flux)/min(flux));
    ANG_BIN(5,kObj)         = 180 - 180/pi*atan(diff(x)/diff(y));
    

    elasped_time = toc(t0)
end

%
path_res = '/home/omartin/Projects/KVS/_results/binaryTests/';
fitswrite(prime_param_piston,[path_res,'Gl569A_prime_param_piston.fits'])
fitswrite(prime_param_piston,[path_res,'Gl569A_prime_param_petal.fits'])
fitswrite(SEP_BIN,[path_res,'Gl569_separationInMas_',num2str(nRes),'pix_ron_',num2str(ron),'.fits'])
fitswrite(DMAG_BIN,[path_res,'Gl569_differentialFlux_',num2str(nRes),'pix_ron_',num2str(ron),'.fits'])
fitswrite(ANG_BIN,[path_res,'Gl569_angle_',num2str(nRes),'pix_ron_',num2str(ron),'.fits'])
fitswrite(im_fit,[path_res,'Gl569_imfit_',num2str(nRes),'pix_ron_',num2str(ron),'.fits'])

%% RESULTS
close all;clc;
suff = '_ron_0';
path_res = '/home/omartin/Projects/KVS/_results/binaryTests/';
nRes        = 100;
SEP_BIN     = fitsread([path_res,'Gl569_separationInMas_100pix_ron_0.fits']);
DMAG_BIN    = fitsread([path_res,'Gl569_differentialFlux_100pix_ron_0.fits']);
ANG_BIN     = fitsread([path_res,'Gl569_angle_100pix_ron_0.fits']);
im_fit      = fitsread([path_res,'Gl569_imfit_100pix_ron_0.fits']);
SEP_BIN_19  = fitsread([path_res,'Gl569_separationInMas_100pix_ron_0_psfAO19.fits']);
DMAG_BIN_19 = fitsread([path_res,'Gl569_differentialFlux_100pix_ron_0_psfAO19.fits']);
ANG_BIN_19  = fitsread([path_res,'Gl569_angle_100pix_ron_0_psfAO19.fits']);
im_fit_19   = fitsread([path_res,'Gl569_imfit_100pix_ron_0_psfAO19.fits']);

figure;
plot(SEP_BIN(1,:),'gs--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','g');
hold on
plot(SEP_BIN(2,:),'mo--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','m');
plot(SEP_BIN(3,:),'cd--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','c');
plot(SEP_BIN_19(:),'rp--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','r');
plot(SEP_BIN(4,:),'bd--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','b');
plot(SEP_BIN(5,:),'k','linewidth',1.2);
legend({'DGC','PSF-R','PETAL-PRIME','PISTON-PSFAO19','PISTON-PRIME','Gl569A'},'interpreter','latex','FontSize',16,'Location','northeast');
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex' );
pbaspect([1,1,1]);
ylabel('Separation (mas)','interpreter','latex','FontSize',20);
xlabel('\# frame','interpreter','latex','FontSize',20);

figure;
plot(DMAG_BIN(1,:),'gs--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','g');
hold on
plot(DMAG_BIN(2,:),'mo--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','m');
plot(DMAG_BIN(3,:),'cd--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','c');
plot(DMAG_BIN_19(:),'rp--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','r');
plot(DMAG_BIN(4,:),'bd--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','b');
plot(DMAG_BIN(5,:),'k','linewidth',1.2);
legend({'DGC','PSF-R','PETAL-PRIME','PISTON-PSFAO19','PISTON-PRIME','Gl569A'},'interpreter','latex','FontSize',16,'Location','northeast');
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex' );
pbaspect([1,1,1]);
ylabel('Differential photometry (mag)','interpreter','latex','FontSize',20);
xlabel('\# frame','interpreter','latex','FontSize',20);

figure;
plot(ANG_BIN(1,:),'gs--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','g');
hold on
plot(ANG_BIN(2,:),'mo--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','m');
plot(ANG_BIN(3,:),'cd--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','c');
plot(ANG_BIN_19(:),'rp--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','r');
plot(ANG_BIN(4,:),'bd--','linewidth',1.2,'markerSize',7,'MarkerFaceColor','b');
plot(ANG_BIN(5,:),'k','linewidth',1.2);
legend({'DGC','PSF-R','PETAL-PRIME','PISTON-PSFAO19','PISTON-PRIME','Gl569A'},'interpreter','latex','FontSize',16,'Location','northeast');
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex' );
pbaspect([1,1,1]);
ylabel('Position angle (degrees)','interpreter','latex','FontSize',20);
xlabel('\# frame','interpreter','latex','FontSize',20);

%% statistics
[median(SEP_BIN(1:3,:),2)',median(SEP_BIN_19)',median(SEP_BIN(4:end,:),2)']
[std(SEP_BIN(1:3,:),[],2)',std(SEP_BIN_19)',std(SEP_BIN(4:end,:),[],2)']

[median(DMAG_BIN(1:3,:),2)',median(DMAG_BIN_19)',median(DMAG_BIN(4:end,:),2)']
[std(DMAG_BIN(1:3,:),[],2)',std(DMAG_BIN_19)',std(DMAG_BIN(4:end,:),[],2)']

[median(ANG_BIN(1:3,:),2)',median(ANG_BIN_19)',median(ANG_BIN(4:end,:),2)']
[std(ANG_BIN(1:3,:),[],2)',std(ANG_BIN_19)',std(ANG_BIN(4:end,:),[],2)']

%%
close all;clc;
k=1;
nBox = 50;
im_fit(:,:,6:7,:) = im_fit(:,:,5:6,:);
im_fit(:,:,5,:)   = im_fit_19;
[yref,xref] = find(im_fit(:,:,1,k) == max(max(im_fit(:,:,1,k))));
idy = yref-nBox/2+1:yref+nBox/2;
idx = yref-nBox/2+1:yref+nBox/2;

A = [im_fit(idy,idx,1,k),im_fit(idy,idx,5,k),im_fit(idy,idx,7,k)];
figure;
imagesc(log10(A.*(A>0)),[1,3.5]);
pbaspect([3,1,1])
set(gca,'XTick',[],'YTick',[]);
cb = colorbar();
cb.FontSize = 20;
cb.TickLabelInterpreter = 'latex';

figure;
subplot(3,1,1)
imagesc(asinh(reshape(im_fit(idy,idx,1,k),nBox,[])),[3,9]);
pbaspect([1,1,1])
set(gca,'XTick',[],'YTick',[]);
cb = colorbar();
cb.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
subplot(3,1,2)
imagesc(asinh(reshape(im_fit(idy,idx,2:4,k),nBox,[])),[3,9]);
set(gca,'XTick',[],'YTick',[]);
cb = colorbar();
cb.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
pbaspect([3,1,1])
subplot(3,1,3)
imagesc(asinh(reshape(im_fit(idy,idx,5:7,k),nBox,[])),[3,9]);
pbaspect([3,1,1])
set(gca,'XTick',[],'YTick',[]);
cb = colorbar();
cb.FontSize = 20;
cb.TickLabelInterpreter = 'latex';

%imagesc(asinh(reshape(im_fit(idy,idx,1:7,k),nBox,[])),[3,9]);

% residual histograms
im_res      = im_fit(:,:,2:7,:) - im_fit(:,:,1,:);
mm          = 1e3;
be          = linspace(-mm,mm,3e2);
res_DGC     = reshape(reshape(squeeze(im_res(idy,idx,1,:)),nBox^2,40),1,[]);
res_psfr    = reshape(reshape(squeeze(im_res(idy,idx,2,:)),nBox^2,40),1,[]);
res_prime_pe= reshape(reshape(squeeze(im_res(idy,idx,3,:)),nBox^2,40),1,[]);
res_psfao19 = reshape(reshape(squeeze(im_res(idy,idx,4,:)),nBox^2,40),1,[]);
res_prime_pis=reshape(reshape(squeeze(im_res(idy,idx,5,:)),nBox^2,40),1,[]);
res_Gl569A  = reshape(reshape(squeeze(im_res(idy,idx,6,:)),nBox^2,40),1,[]);

figure;
histogram(res_DGC,'Normalization','Probability','BinEdges',be);hold on;
histogram(res_psfr,'Normalization','Probability','BinEdges',be);
histogram(res_prime_pe,'Normalization','Probability','BinEdges',be);
histogram(res_psfao19,'Normalization','Probability','BinEdges',be);
histogram(res_prime_pis,'Normalization','Probability','BinEdges',be);
histogram(res_Gl569A,'Normalization','Probability','BinEdges',be);
legend({'DGC','PSF-R','PETAL-PRIME','PISTON-PSFAO19','PISTON-PRIME','Gl569A'},'interpreter','latex','FontSize',16,'Location','northeast');
%histfit(res_DGC)
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex' );
pbaspect([1,1,1]);
ylabel('Probability','interpreter','latex','FontSize',20);
xlabel('Residual intensity (ADU)','interpreter','latex','FontSize',20);
xlim([-mm,mm])

% fit distributions with a gaussian function
pd = fitdist(res_DGC','Normal');
sig_dist(1) = pd.sigma;
med_dist(1) = pd.median;
pd = fitdist(res_psfr','Normal');
sig_dist(2) = pd.sigma;
med_dist(2) = pd.median;
pd = fitdist(res_prime_pe','Normal');
sig_dist(3) = pd.sigma;
med_dist(3) = pd.median;
pd = fitdist(res_psfao19','Normal');
sig_dist(4) = pd.sigma;
med_dist(4) = pd.median;
pd = fitdist(res_prime_pis','Normal');
sig_dist(5) = pd.sigma;
med_dist(5) = pd.median;
pd = fitdist(res_Gl569A','Normal');
sig_dist(6) = pd.sigma;
med_dist(6) = pd.median;

%% OLD
   
    %8\ Measure photometry/astrometry with PRIME
%     p.getPrimePSF('objname',objname{kObj},'resolution',nRes,'fov',nRes+20,'ron',300,...
%         'fitStatModes',1:36,'statModesFunction','piston','fitGains',[true,true,false],'umax',umax);
%     y                       = p.psfp.catalog_fit.y;
%     x                       = p.psfp.catalog_fit.x;
%     flux                    = p.psfp.catalog_fit.flux;
%     SEP_BIN(5,kObj)         = hypot(diff(x),diff(y));
%     DMAG_BIN(5,kObj)        = 2.5*log10(max(flux)/min(flux));
%     ANG_BIN(5,kObj)         = 180 - 180/pi*atan(diff(y)/diff(x));

    
    % prime PSF with anisoplanatism
%     if includeAniso
%         p.trs.cam.resolution    = 150;
%         xdata{1} = 1:p.trs.atm.nLayer;
%         xdata{2} = 1:p.trs.atm.nLayer;
%         xdata{3} = p.trs.atm.nLayer+1;
%         xdata{4} = p.trs.atm.nLayer+2;
%         xdata{6} = p.trs.atm.nLayer+3:p.trs.atm.nLayer+2+36;
%         p.trs.ngs.x         = 0;
%         p.trs.ngs.y         = 0;
%         p.trs.src           = p.trs.src(1);
%         p.psfr              = computeAnisoplanatismPhaseStructureFunction(p.psfr);
%         p.psfp.weightMap    = 1;
%         p.trs.src(iBin).nObj = 1;
%         Cn2                 = p.psfp.x_fixed.Cn2;
%         psf_prime           = imageModel([1,0,0,Cn2*p.psfp.x_final(4)/sum(Cn2),prime_param_piston(5:end-1,kObj)',0],xdata,p.psfp);
%         psf_prime           = psf_prime/sum(psf_prime(:));
%         p.trs.cam.resolution    = nRes;
%     end
   