
function Dani_l = instantiateAnisoplanatism(psfr,gs,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfReconstruction'));
inputs.addRequired('gs',@isstruct);
inputs.addParameter('isTT',false,@islogical);
inputs.parse(psfr,gs,varargin{:});


fprintf('Instantiate the anisoplanatism model...\n');

%1\ Defining the spatial filters
D    = psfr.trs.tel.Ddm;
npt  = psfr.otf.dk;
isTT = inputs.Results.isTT;
atm  = psfr.trs.atm;
nSrc = numel(psfr.trs.src);

if isTT % -> anisokinetism as a Gaussian kernel

    % defining tip-tilt modes 
    zern = zernike(2:3,'resolution',psfr.otf.nOtf);
    X    = reshape(zern.modes(:,1),psfr.otf.nOtf,psfr.otf.nOtf);
    Y    = reshape(zern.modes(:,2),psfr.otf.nOtf,psfr.otf.nOtf);
    X2   = X.^2;
    Y2   = Y.^2;
    XY   = X.*Y';
    YX   = Y.*X';
    % defining the oomao-like atmosphere
    atm_oomao = atmosphere(atm.wavelength,atm.r0,atm.L0,'fractionnalR0',atm.weights,'altitude',atm.heights,'windSpeed',atm.windSpeed,'windDirection',atm.windDirection);
    % defining oomao-like sources
    [aS,zS]    = cart2pol(gs.x,gs.y);
    ngs        = source('zenith',zS,'azimuth',aS);
    [aS,zS]    = cart2pol(psfr.trs.src.x,psfr.trs.src.y);
    sci        = source('zenith',zS,'azimuth',aS); 
    % instantiating the phase structure function
    Dani_l = zeros(psfr.otf.nOtf,psfr.otf.nOtf,atm.nLayer,nSrc);
    
    % computing the phase structure function for each layers and src
    for iSrc = 1:nSrc
        thx = (psfr.trs.src(iSrc).x - gs.x);
        thy = (psfr.trs.src(iSrc).y - gs.y);
        if thx ~=0 || thy~=0
            for l = 1:atm.nLayer
                if atm.heights(l) ~=0
                    % update the atmosphere
                    atm_l = atm_oomao.slab(l);
                    fr0   = atm_l.r0^(-5/3);% * atm_oomao.layer(l).fractionnalR0;
                    % Get the 2x2 covariance matrices of the tip-tilt
                    covTT = zernikeStats.angularCovariance(zern,atm_l,[sci(iSrc),sci(iSrc),ngs]);
                    % get the 2x2 anisokinetism covariance matrix
                    covAniso = covTT{1,1} + covTT{3,3} - covTT{1,3} - covTT{3,1};
                    % defining the Gaussian kernel
                    Dani_l(:,:,l,iSrc) = (covAniso(1,1) * X2 + covAniso(2,2)*Y2 + covAniso(1,2)*XY + covAniso(2,1)*YX)/fr0;
                end
            end
        end
    end
    %x      = (-1+1/npt:2/npt:1-1/npt);
    %[X,Y]  = meshgrid(x,x);
    %TT     = [X(:),Y(:)];
    %zern   = zernike_puako(2:3,npt);
    %TT     = zern.modes;
    %Hfilter= TT*pinv(TT);
    
else %-> focal-angular anisokinetism

    if gs.height ~= Inf % focal anisoplanatism -> TT filtering
        zern   = zernike_puako(2:3,npt);
        TT     = zern.modes;
        Hfilter= eye(npt^2) - TT*pinv(TT);
        Hfilter= Hfilter*psfr.trs.mat.Hdm;
    else
        Hfilter= psfr.trs.mat.Hdm;
    end
    
    %2\ SF Calculation
    if psfr.flags.toeplitz
        Dani_l = zeros(psfr.otf.nOtf,psfr.otf.nOtf,atm.nLayer,nSrc);
    else
        Dani_l = zeros(psfr.otf.dk^2,psfr.otf.dk^2,atm.nLayer,nSrc);
    end

    if strcmpi(psfr.flags.anisoMethod,'FLICKER') % Use the Flicker's 2008 report to derive the SF
        % Get inputs
        f0      = 2*pi/atm.L0;
        % Phase sample locations in the pupil
        x       = -D/2:D/(psfr.otf.dk-1):D/2;
        [x1,y1] = meshgrid(x);
        X1      = (x1(:)*ones(1,psfr.otf.dk^2))';
        Y1      = repmat(y1,[psfr.otf.dk,psfr.otf.dk]);
        % Samples separation in the pupil
        rhoX    = bsxfun(@minus,X1,x1(:));
        rhoY    = bsxfun(@minus,Y1,y1(:)');
        % Instantiation
        Ialpha  = @(x,y) mcDonald(f0*hypot(x,y));
        I0      = Ialpha(0,0);
        I1      = Ialpha(rhoX,rhoY);
        cte     = 0.12184*0.06*(2*pi)^2;
        
        % Anisoplanatism Structure Function
        %rad2arcsec = 3600 * 180/pi;
        for iSrc = 1:nSrc
            thx = (psfr.trs.src(iSrc).x - gs.x);
            thy = (psfr.trs.src(iSrc).y - gs.y);
            
            for l = 1:atm.nLayer
                zl   = atm.heights(l);
                if zl > 0
                    if gs.height == Inf  || isempty(gs.height)
                        I2    = Ialpha(rhoX+zl*thx,rhoY+zl*thy);
                        I3    = Ialpha(zl*thx,zl*thy);
                        I4    = Ialpha(rhoX-zl*thx,rhoY-zl*thy);
                        tmp   = Hfilter*(2*I0 - 2*I1 + I2 - 2*I3  + I4)*Hfilter';
                    else
                        g     = zl/gs.height;
                        I2    = Ialpha(rhoX*(1-g) , rhoY*(1-g));
                        I3    = Ialpha(rhoX -g*X1 + zl*thx , rhoY - g*Y1 + zl*thy);
                        I4    = Ialpha(g*X1 - zl*thx , g*Y1 - zl*thy);
                        I5    = Ialpha(g*(rhoX-X1) -zl*thx , g*(rhoY-Y1) - zl*thy);
                        I6    = Ialpha((1-g)*rhoX + g*X1 - zl*thx , (1-g)*rhoY + g*Y1 - zl*thy);
                        tmp   = Hfilter*(2*I0 - I1 - I2 + I3 - I4 - I5 + I6)*Hfilter';
                    end
                    
                    if psfr.flags.toeplitz
                        tmp = puakoTools.covMatrix2Map(tmp,psfr.otf.dk,psfr.otf.dk);
                        tmp = puakoTools.interpolateOtf(tmp,psfr.otf.nOtf);
                    end
                    Dani_l(:,:,l,iSrc)  = Dani_l(:,:,l,iSrc) + cte*atm.L0^(5/3)*tmp;
                end
            end
        end
        
    else % Use the native OOMAO routines in phaseStats to calculate the anisoplanatism SF.
        
        %Layers-by-layers covariance matrices
        cov             = @(dk,D,atm,gs) puakoTools.spatioAngularCovarianceMatrix(dk,D,atm,gs);
        cross           = @(dk,D,atm,src,gs) puakoTools.spatioAngularCovarianceMatrix(dk,D,atm,src,'srcCC',gs);
        thereIsFocAniso = gs.height~=Inf;
        
        for l=1:atm.nLayer
            if atm.heights(l) ~=0 % Check that the layer is not at the ground
                subAtm.r0               = 1;
                subAtm.L0               = atm.L0;
                subAtm.heights          = atm.heights(l);
                subAtm.weights          = 1;
                subAtm.windSpeed        = atm.windSpeed(l);
                subAtm.windDirection    = atm.windDirection(l);
                subAtm.nLayer           = 1;
                
                % Atmospheric covariance matrix in the guide star direction
                Cgg = cov(psfr.otf.dk,D,subAtm,gs);
                
                for iSrc=1:psfr.trs.src.nSrc
                    thereIsAngAniso = any([psfr.trs.src(iSrc).x - gs.x,psfr.trs.src(iSrc).y - gs.y]);
                    % Phase covariance matrix
                    if psfr.trs.src(iSrc).height ~= gs(iSrc).height
                        %LGS case: Css ~= Csg
                        Css = cov(psfr.otf.dk,D,subAtm,psfr.trs.src(iSrc));
                    else
                        %NGS case: Css = Csg
                        Css = Cgg;
                    end
                    % Cross covariance matrix
                    if thereIsAngAniso || thereIsFocAniso
                        Cgs = cross(psfr.otf.dk,D,subAtm,psfr.trs.src(iSrc),gs);
                        Cgs = Cgs{1};
                    else
                        Cgs = 0;
                    end
                    % Anisoplanatic covariance matrix
                    Cani = Hfilter*(Css + Cgg - Cgs - Cgs')*Hfilter';
                    if psfr.flags.toeplitz
                        Cani = puakoTools.covMatrix2Map(Cani,psfr.otf.dk,psfr.otf.dk);
                        Cani = interpolateOtf(Cani,psfr.otf.nOtf);
                        Dani_l(:,:,l,iSrc) = 2*(max(Cani(:)) - Cani);
                    else
                        Dani_l(:,:,l,iSrc) = 2*(diag(Cani) - Cani);
                    end
                end
            end
        end
    end
end

function out = mcDonald(x)
        out = x.^(5/6.).*besselk(5./6,x)./(2^(5/6.)*gamma(11/6.)) ;
        out(find(x == 0)) = 3/5.;
        
function out = interpolateOtf(otf,nRes,method)
            if nargin < 3
                method = 'spline';
            end
            
            if numel(nRes) == 1
                nRes = [nRes nRes];
            end
            
            % Define input meshgrid
            notf_x = size(otf,1);
            notf_y = size(otf,2);
            if mod(notf_x,2)==0
                u1D = (-notf_x/2:1:notf_x/2-1)*2/notf_x;
            else
                u1D = (-floor(notf_x/2):1:floor(notf_x/2))*2/notf_x;
            end
            if mod(notf_y,2)==0
                v1D = (-notf_y/2:1:notf_y/2-1)*2/notf_y;
            else
                v1D = (-floor(notf_y/2):1:floor(notf_y/2))*2/notf_y;
            end
            [Uxi,Uyi]= meshgrid(u1D,v1D);
            
            % Define output meshgrid
            nRes_x = nRes(1);
            nRes_y = nRes(2);
            
            if mod(nRes_x,2)==0
                u1D2 = (-nRes_x/2:1:nRes_x/2-1)*2/nRes_x;
            else
                u1D2 = (-floor(nRes_x/2):1:floor(nRes_x/2))*2/nRes_x;
            end
            
            if mod(nRes_y,2)==0
                v1D2 = (-nRes_y/2:1:nRes_y/2-1)*2/nRes_y;
            else
                v1D2 = (-floor(nRes_y/2):1:floor(nRes_y/2))*2/nRes_y;
            end
            [Uxo,Uyo]= meshgrid(u1D2,v1D2);
            
            % Interpolation
            out = interp2(Uxi,Uyi,otf,Uxo,Uyo,method);
        

