function displayResults(obj,varargin)
inputs = inputParser;
inputs.addRequired('obj', @(x) isa(x,'psfReconstruction') || isa(x,'PRIME')  );
inputs.addParameter('fontsize', 20,@isnumeric );
inputs.parse(obj,varargin{:});

fontsize = inputs.Results.fontsize;

if isa(obj,'psfReconstruction')
    resultsof = 'PSF-R';
    ps = obj(1).trs.cam.pixelScale;
    samp = obj(1).trs.cam.samp;
else
    resultsof = 'PRIME';
    ps = obj.psfr.trs.cam.pixelScale;
    samp = obj.psfr.trs.cam.samp;
end




if numel(obj) == 1 % Results vizualization of a single case
  
    %1\ Grab images
    im_sky = obj.sky.image;
    im_rec = obj.psf.im_fit; % the reconstructed PSF fitted over the sky to get same scale and astrometry
    im_sky(im_sky<0)=0;
    im_res = im_rec - im_sky;
    
    fov  =size(im_sky,1);
    rov  = floor(fov/2);
    [~,~,ron] = tools.getFlux(im_sky);
    mx = max(max(im_sky(:)),max(im_rec(:)));
    mn = min([ron,min(radial(im_rec,rov+1,rov+1)),min(radial(im_sky,rov+1,rov+1))]);
      
    %2\ 1D plot of PSF
    x_r = linspace(0,ps*fov/2,rov);
    x_f = ps*rov*linspace(-1,1,fov);

    
    h = figure;
    semilogy(x_r,radial(im_sky,rov+1,rov+1),'b-');hold on;
    semilogy(x_r,radial(im_rec,rov+1,rov+1),'r--');
    semilogy(x_r,abs(radial(im_res,rov+1,rov+1)),'k:');
    ylim([mn,mx]);
    xlabel('Angular distance (mas)','interpreter','latex','FontSize',fontsize);
    ylabel('Azimuthal profile (ADU)','interpreter','latex','FontSize',fontsize);
    legend({'Sky image',resultsof,'Residual'},'interpreter','latex','FontSize',fontsize,'Location','northeast');
    set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex' );
    pbaspect([1.6,1,1]);
    
    h = figure;
    semilogy(x_f,im_sky(:,rov+1),'b-');hold on;
    semilogy(x_f,im_rec(:,rov+1),'r--');
    semilogy(x_f,abs(im_res(:,rov+1)),'k:');
    ylim([mn,mx]);
    xlabel('Angular distance (mas)','interpreter','latex','FontSize',fontsize);
    ylabel('x-axis profile (ADU)','interpreter','latex','FontSize',fontsize);
    legend({'Sky image',resultsof,'Residual'},'interpreter','latex','FontSize',fontsize,'Location','northeast');
    set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex' );
    pbaspect([1.6,1,1]);
    
    h = figure;
    semilogy(x_f,im_sky(rov+1,:),'b-');hold on;
    semilogy(x_f,im_rec(rov+1,:),'r--');
    semilogy(x_f,abs(im_res(rov+1,:)),'k:');
    ylim([mn,mx]);
    legend({'Sky image',resultsof,'Residual'},'interpreter','latex','FontSize',fontsize,'Location','northeast');
    xlabel('Angular distance (mas)','interpreter','latex','FontSize',fontsize);
    ylabel('y-axis profile (ADU)','interpreter','latex','FontSize',fontsize);
    set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex' );
    pbaspect([1.6,1,1]);

    %2\ 1D plot of OTF
    otf_sky = abs(fftshift(fft2(im_sky)));
    otf_sky = otf_sky/max(otf_sky(:));
    otf_rec = abs(fftshift(fft2(im_rec)));
    otf_rec = otf_rec/max(otf_rec(:));
    otf_res = otf_rec - otf_sky;
    
    u       = linspace(0,max(samp,1),fov/2);
    h = figure;
    semilogy(u,(radial(otf_sky)),'b--');hold on;
    semilogy(u,(radial(otf_rec)),'r-');
    semilogy(u,abs(radial(otf_res)),'k:');
    xlabel('Angular frequency (D/$\lambda$)','interpreter','latex','FontSize',fontsize);
    ylabel('Modulation transfer function amplitude','interpreter','latex','FontSize',fontsize);
    xlim([0,1]);
    legend({'Sky MTF',resultsof,'Residual'},'interpreter','latex','FontSize',20);
    set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex' );
    pbaspect([1.6,1,1]);

    %3\ 2D PSF;
    P = [im_sky,im_rec,im_res];
    Plog = log10(abs(P));
    if sign(max(Plog(:))) ==1
        mn = 0.2*max(Plog(:));
        mx = 0.9*max(Plog(:));
    else
        mn = 5*max(Plog(:));
        mx = 1.1*min(Plog(:));
    end
    
    h = figure;
    imagesc(Plog,[mn,mx]);
    set(gca,'XTick',[],'YTick',[]);
    pbaspect([3,1,1])
else
    
    sky = [obj.sky];
    rec = [obj.psf];
    
    %1\ Compare  Strehl values
    SRsky = 1e2*[sky.SR];
    dSRsky = 1e2*[sky.dSR];
    SR = 1e2*[rec.SR];
        
    h = figure;
    errorbar(SRsky,SR,[],[],-dSRsky/2,dSRsky/2,'bs','MarkerFaceColor','b','MarkerSize',5);
    hold on;
    xx = [xlim()];
    yy = [ylim()];
    mn = min(min([xx(:),yy(:)]));
    mx = max(max([xx(:),yy(:)]));
    plot([mn,mx],[mn,mx],'k--');
    ylabel('Reconstructed Strehl ratio (\%)','interpreter','latex','FontSize',fontsize);
    xlabel('Image Strehl ratio (\%)','interpreter','latex','FontSize',fontsize);
    set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex');
    pbaspect([1.6,1,1]);
    
     %2\ Compare FWHM values
     
    FWHMxsky = [sky.FWHMx];
    FWHMysky = [sky.FWHMy];
    dFWHMsky = [sky.dFWHM];
    FWHMx = [rec.FWHMx];
    FWHMy = [rec.FWHMy];
    
    h = figure;
    errorbar(FWHMxsky,FWHMx,-dFWHMsky/2,dFWHMsky/2,-dFWHMsky/2,dFWHMsky/2,'bs','MarkerFaceColor','b','MarkerSize',5);
    hold on;
    errorbar(FWHMysky,FWHMy,-dFWHMsky/2,dFWHMsky/2,-dFWHMsky/2,dFWHMsky/2,'rs','MarkerFaceColor','r','MarkerSize',5);
     xx = [xlim()];
    yy = [ylim()];
    mn = min(min([xx(:),yy(:)]));
    mx = max(max([xx(:),yy(:)]));
    plot([mn,mx],[mn,mx],'k--');
    ylabel('Reconstructed FWHM (%)','interpreter','latex','FontSize',fontsize);
    xlabel('Image FWHM (%)','interpreter','latex','FontSize',fontsize);
    legend({'X-axis','Y-axis'},'interpreter','latex','Fontsize',fontsize,'Location','northeast');
    set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex');
    pbaspect([1.6,1,1]);


end





