function displayResults(obj,varargin)
inputs = inputParser;
inputs.addRequired('obj', @(x) isa(x,'psfReconstruction') || isa(x,'prime')  );
inputs.addParameter('fov', size(obj.psf.image,1),@isnumeric );
inputs.addParameter('fontsize', 16,@isnumeric );
inputs.addParameter('linewidth', 1.5,@isnumeric );
inputs.parse(obj,varargin{:});
fov         = inputs.Results.fov;
fontsize    = inputs.Results.fontsize;
linewidth   = inputs.Results.linewidth;


if isa(obj,'psfReconstruction')
    resultsof = 'PSF-R';
    ps      = obj(1).trs.cam.pixelScale;
    samp    = obj(1).trs.cam.samp;
    trs     = obj.trs; 
    im_rec  = obj.psf(1).im_fit; % the reconstructed PSF fitted over the sky to get same scale and astrometry
    [x0,y0] = puakoTools.cog(im_rec);
    x0      = round(x0);y0 = round(y0);
else
    resultsof = 'PRIME';
    ps      = obj.psfr.trs.cam.pixelScale;
    samp    = obj.psfr.trs.cam.samp;
    trs     = obj.psfr.trs;
    im_rec  = obj.psf(1).image;
    x0      = obj.catalog_fit.x(1)/trs.cam.pixelScale + fov/2;
    y0      = obj.catalog_fit.y(1)/trs.cam.pixelScale + fov/2;
end


%1\ Grab images
im_sky = trs.sky.image;
if fov~=size(obj.psf.image,1)
    im_sky = puakoTools.crop(im_sky,fov);
    im_rec = puakoTools.crop(im_rec,fov);
end

F = max(im_sky(:));
im_sky = im_sky/F;
im_rec = im_rec/F;
im_res = abs(im_rec - im_sky);

fov  =size(im_sky,1);
rov  = floor(fov/2);
[~,~,ron] = puakoTools.getFlux(im_sky);
mx = 1.1*max(max(im_sky(:)),max(im_rec(:)));
mn = min([ron,min(radial(im_rec,rov+1,rov+1)),min(radial(im_sky,rov+1,rov+1))]);

%2\ 1D plot of PSF
x_r = linspace(0,ps*fov/2,rov);
x_f = ps*rov*linspace(-1,1,fov);

h = figure;
subplot(2,2,1)
plot(x_r,log10(radial(im_sky,x0,y0)),'b-','linewidth',linewidth);hold on;
plot(x_r,log10(radial(im_rec,x0,y0)),'r','linewidth',linewidth);
plot(x_r,log10(radial(im_res,x0,y0)),'k:','linewidth',linewidth);
ylim([-5,0.5]);
xlim([min(x_r)-10,max(x_r(:))])
xlabel('Angular distance (mas)','interpreter','latex','FontSize',fontsize);
ylabel('Azimuthal profile (ADU)','interpreter','latex','FontSize',fontsize);
legend({'Sky image',resultsof,'Residual'},'interpreter','latex','FontSize',fontsize,'Location','northeast');
set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex' );
pbaspect([1.6,1,1]);

subplot(2,2,3)
plot(x_f,log10(im_sky(:,round(y0))),'b-','linewidth',linewidth);hold on;
plot(x_f,log10(im_rec(:,round(y0))),'r','linewidth',linewidth);
plot(x_f,log10(abs(im_res(:,round(y0)))),'k:','linewidth',linewidth);
ylim([-5,0.5]);
xlim([min(x_f),max(x_f(:))])
xlabel('Angular distance (mas)','interpreter','latex','FontSize',fontsize);
ylabel('x-axis profile (ADU)','interpreter','latex','FontSize',fontsize);
legend({'Sky image',resultsof,'Residual'},'interpreter','latex','FontSize',fontsize,'Location','northeast');
set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex' );
pbaspect([1.6,1,1]);

subplot(2,2,4)
plot(x_f,log10(im_sky(round(x0),:)),'b-','linewidth',linewidth);hold on;
plot(x_f,log10(im_rec(round(x0),:)),'r','linewidth',linewidth);
plot(x_f,log10(im_res(round(x0),:)),'k:','linewidth',linewidth);
ylim([-5,0.5]);
xlim(([min(x_f),max(x_f(:))]))
legend({'Sky image',resultsof,'Residual'},'interpreter','latex','FontSize',fontsize,'Location','northeast');
xlabel('Angular distance (mas)','interpreter','latex','FontSize',fontsize);
ylabel('y-axis profile (ADU)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex' );
pbaspect([1.6,1,1]);

%2\ 1D plot of OTF
subplot(2,2,2)
otf_rec = puakoTools.psf2otf(im_rec);
mtf_rec = abs(otf_rec/max(otf_rec(:)));
otf_sky = puakoTools.psf2otf(im_sky);
mtf_sky = abs(otf_sky/max(otf_rec(:)));
mtf_res = abs(mtf_rec - mtf_sky);
mtf_tel = puakoTools.pupil2otf(trs.tel.pupil,trs.tel.pupil*0,1);
mtf_tel = puakoTools.interpolateOtf(mtf_tel,(fov/trs.cam.samp));
mtf_tel = puakoTools.enlargePupil(mtf_tel,trs.cam.samp);
mtf_tel = puakoTools.interpolateOtf(mtf_tel,(fov));

u       = linspace(0,max(samp,1),(fov/2));
[ym,xm] = find(mtf_tel == max(mtf_tel(:)));
mtel    = (radial(mtf_tel,ym,xm));
semilogy(u,mtel,'k--','linewidth',linewidth);hold on;
ms      = (radial(mtf_sky,fov/2+1,fov/2+1));
semilogy(u,ms,'b-','linewidth',linewidth);
mr      = (radial(mtf_rec,fov/2+1,fov/2+1));
semilogy(u,mr,'r-','linewidth',linewidth);
dr      = radial(mtf_res,fov/2+1,fov/2+1);
semilogy(u,dr,'k:','linewidth',linewidth);
xlabel('Angular frequency (D/$\lambda$)','interpreter','latex','FontSize',fontsize);
ylabel('MTF profile','interpreter','latex','FontSize',fontsize);
xlim([0,1]);
id = find(u>1);
ylim([1e-5,1.05])
legend({'Telescope MTF','Sky MTF',resultsof,'Residual'},'interpreter','latex','FontSize',fontsize,'Location','southwest');
set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex' );
pbaspect([1.6,1,1]);

%3\ 2D PSF;
P = [im_sky,im_rec,im_res];
Plog = log10(P);
h = figure;
imagesc(Plog,[-5,0.5]);
set(gca,'XTick',[],'YTick',[]);
pbaspect([3,1,1])
cb = colorbar();
cb.TickLabelInterpreter = 'latex';
cb.FontSize = fontsize;

%4\ Static maps;
if strcmp(resultsof,'PRIME') && ~isempty(obj.map_fit) && ~isempty(obj.map_fit.coefs)
    f = obj.psfr.trs.cam.wavelength*1e9/2/pi;
    figure;
    imagesc(f*[obj.map_fit.map_stat-obj.map_fit.map_zer,obj.map_fit.map_zer,obj.map_fit.map_stat])
    set(gca,'XTick',[],'YTick',[]);
    pbaspect([3,1,1])
    cb = colorbar();
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = fontsize;
    
    figure;
    plot(obj.map_fit.jindex,obj.map_fit.coefs,'bs--','MarkerFaceColor','b','MarkerSize',7);
    xlabel('Mode index','interpreter','latex','fontsize',fontsize);
    ylabel('Modal coefficients (nm)','interpreter','latex','fontsize',fontsize);
    set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex' );
end

end





