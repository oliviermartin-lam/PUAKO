%{
------------HEADER-----------------
Objective          ::  Compute bi-dimensional phase structure function of the anisoplanatism error

INPUT VARS
 psfr          :: The psfReconstruction top-level class
 
OUTPUT VARS
 Dani_l             :: bi-dimensional layer-by-layer phase structure function map  (psfr.nOtf x psfr.nOtf x psfr.atm.nLayer x psfr.nStars)
Dani                  :: bi-dimensional integrated phase structure function map  (psfr.nOtf x psfr.nOtf x psfr.nStars)
Kani                  :: Anisoplanatism spatial filter based on psfr.atm (psfr.nOtf x psfr.nOtf x psfr.nStars)
Created by         :: O. Beltramo-Martin - ONERA/LAM
Creation date      :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function [Dani_l,Dani,Kani] = computeAnisoplanatismPhaseStructureFunction(psfr)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfReconstruction'));
inputs.parse(psfr);


%1\ Verify if there is anisoplanatism
Dani_l = 0;Dani = 0;Kani = 1;

psfr.flags.isFocalAniso = strcmp(psfr.trs.aoMode,'LGS');
psfr.flags.isAngularAniso = strcmp(psfr.trs.aoMode,'NGS') && (psfr.trs.ngs.x~=psfr.trs.src.x || psfr.trs.ngs.y~=psfr.trs.src.y);
psfr.flags.isTiltAniso = strcmp(psfr.trs.aoMode,'LGS') && (psfr.trs.ngs.x~=0 || psfr.trs.ngs.y~=0);    
psfr.flags.isAniso = psfr.flags.isFocalAniso + psfr.flags.isAngularAniso + psfr.flags.isTiltAniso;

if psfr.flags.isAniso    
    Cn2 = psfr.trs.atm.Cn2/sum(psfr.trs.atm.Cn2(:))*psfr.trs.res.seeing.r0*(psfr.trs.cam.wavelength/0.5e-6)^1.2;
    %2\ Get the anisoplanatism phase structure function due to angular and focal anisoplanatism
    if psfr.flags.isAngularAniso
    psfr.sf.Dani_l= instantiateAnisoplanatism(psfr,psfr.trs.ngs);   
    elseif psfr.flags.isFocalAniso
        psfr.sf.Dani_l= instantiateAnisoplanatism(psfr,psfr.trs.lgs);
    end
    psfr.sf.Dani  = squeeze(sum(bsxfun(@times, psfr.sf.Dani_l, reshape(Cn2,1,1,[])), 3));
    
    if psfr.flags.isTiltAniso
    %3\ Get the tip-tilt anisoplanatism phase structure function
        psfr.sf.DaniTT_l  = instantiateAnisoplanatism(psfr,psfr.trs.ngs,'isTT',true);
        psfr.sf.Dani  = squeeze(psfr.sf.Dani+ sum(bsxfun(@times, psfr.sf.DaniTT_l, reshape(Cn2,1,1,[])), 3));
        psfr.sf.Dani_l= psfr.sf.Dani_l  + psfr.sf.DaniTT_l;
    end
    
    %4\ Get the anisoplanatism spatial filter
    if psfr.flags.toeplitz
        Kani = exp(-0.5*psfr.sf.Dani);
    else
        idx       = true(psfr.otf.dk);
        dp        = psfr.tel.Ddm/(psfr.otf.dk-1);
        psfr.otfDen= tools.zonalCovarianceToOtf(Dani*0,psfr.otf.nOtf,psfr.trs.tel.Ddm,dp,idx);
        psfr.mskOtf= psfr.otfDen > 1e-6;
        Kani  = tools.zonalCovarianceToOtf(Dani,psfr.otf.nOtf,psfr.trs.tel.Ddm,dp,idx);
        Kani(psfr.mskOtf) = Kani(psfr.mskOtf)./psfr.otfDen(psfr.mskOtf);
    end
end

function Dani_l = instantiateAnisoplanatism(psfr,gs,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfReconstruction'));
inputs.addRequired('gs',@isstruct);
inputs.addParameter('isTT',false,@islogical);
inputs.parse(psfr,gs,varargin{:});


fprintf('Instantiate the anisoplanatism model...');

%1\ Defining the spatial filters
D       = psfr.trs.tel.Ddm;
npt  = psfr.otf.dk;
isTT = inputs.Results.isTT;
atm = psfr.trs.atm;
if gs.height ~= Inf
    zern   = zernike(2:3,npt);
    TT     = zern.modes;
    Hfilter= eye(npt^2) - TT*pinv(TT);
    Hfilter= Hfilter*psfr.trs.mat.Hdm;
elseif gs.height == Inf && isTT
    x      = (-1+1/npt:2/npt:1-1/npt);
    [X,Y]  = meshgrid(x,x);
    TT     = [X(:),Y(:)];
    Hfilter= TT*pinv(TT);
else
    Hfilter      = psfr.Hdm;
end

%2\ SF Calculation
if psfr.flags.toeplitz
    Dani_l = zeros(psfr.otf.nOtf,psfr.otf.nOtf,atm.nLayer,psfr.trs.src.nSrc);
else
    Dani_l  = zeros(psfr.otf.dk^2,psfr.otf.dk^2,atm.nLayer,psfr.trs.src.nSrc);
end

if strcmpi(psfr.flags.anisoMethod,'FLICKER') % Use the Flicker's 2008 report to derive the SF
    % Get inputs
    f0      = 2*pi/atm.L0;
    % Phase sample locations in the pupil
    x       = -D/2:D/(psfr.dk-1):D/2;
    [x1,y1] = meshgrid(x);
    X1      = (x1(:)*ones(1,psfr.dk^2))';
    Y1      = repmat(y1,[psfr.dk,psfr.dk]);
    % Samples separation in the pupil
    rhoX    = bsxfun(@minus,X1,x1(:));
    rhoY    = bsxfun(@minus,Y1,y1(:)');
    % Instantiation
    Ialpha  = @(x,y) mcDonald(f0*hypot(x,y));
    I0      = Ialpha(0,0);
    I1      = Ialpha(rhoX,rhoY);
    cte     = 0.12184*0.06*(2*pi)^2;
    
    % Anisoplanatism Structure Function
    for iSrc = 1:numel(psfr.trs.src.nSrc)        
        thx = psfr.trs.src(iSrc).x - gs.x;
        thy =psfr.trs.src(iSrc).y - gs.y;
        
        for l = 1:atm.nLayer
            zl   = atm.heights(l);
            if gs.height == Inf  || isempty(gs.height)
                I2    = Ialpha(rhoX+zl*thx,rhoY+zl*thy);
                I3    = Ialpha(zl*thx,zl*thy);
                I4    = Ialpha(rhoX-zl*thx,rhoY-zl*thy);
                tmp   = I0 - 2*I1 + I2 - 2*I3  + I4;
            else
                g     = zl/gs.height;
                I2    = Ialpha(rhoX*(1-g),rhoY*(1-g));
                I3    = Ialpha(rhoX-g*X1+zl*thx,rhoY-g*Y1+zl*thy);
                I4    = Ialpha(g*X1-zl*thx,g*Y1-zl*thy);
                I5    = Ialpha(g*(rhoX-X1)-zl*thx,g*(rhoY-Y1)-zl*thy);
                I6    = Ialpha((1-g)*rhoX+g*X1-zl*thx,(1-g)*rhoY+g*Y1-zl*thy);
                tmp   = 2.*I0 - I1 - I2 + I3 - I4 - I5 + I6;
            end
            
            if psfr.flagToeplitz
                tmp = tools.covMatrix2Map(tmp,psfr.otf.dk,psfr.otf.dk);
                tmp = tools.interpolateOtf(tmp,psfr.otf.nOtf);
            end
            Dani_l(:,:,l,iSrc)  = Dani_l(:,:,l,iSrc) + cte*atm.L0^(5/3)*Hfilter*tmp*Hfilter';
        end
    end
    
else % Use the native OOMAO routines in phaseStats to calculate the anisoplanatism SF.
    
    %Layers-by-layers covariance matrices
    atm.r0 = psfr.trs.res.seeing.r0*(psfr.trs.cam.wavelength/0.5e-6)^1.2;
    cov = @(dk,D,atm,gs) tools.spatioAngularCovarianceMatrix(dk,D,atm,gs);
    cross = @(dk,D,atm,src,gs) tools.spatioAngularCovarianceMatrix(dk,D,atm,src,'srcCC',gs);        
    thereIsFocAniso = gs.height~=Inf;
    Cn2 = psfr.trs.atm.Cn2/sum(psfr.trs.atm.Cn2(:))*psfr.trs.atm.r0^(-5/3);
    
    for l=1:atm.nLayer
        if atm.heights(l) ~=0 % Check that the layer is not at the ground
            f0 =Cn2(l);
            subAtm.r0 = atm.r0^(-5/3)*atm.weights(l);
            subAtm.L0 = atm.L0;
            subAtm.heights = atm.heights(l);
            subAtm.weights = 1;
            subAtm.windSpeed = atm.windSpeed(l);
            subAtm.windDirection = atm.windDirection(l);
            subAtm.nLayer = 1;
            
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
                Cani = Hfilter*(Css + Cgg - Cgs - Cgs')*Hfilter'/f0;
                if psfr.flags.toeplitz
                    Cani = tools.covMatrix2Map(Cani,psfr.otf.dk,psfr.otf.dk);
                    Cani = tools.interpolateOtf(Cani,psfr.otf.nOtf);
                    Dani_l(:,:,l,iSrc) = 2*(max(Cani(:)) - Cani);
                else                
                    Dani_l(:,:,l,iSrc) = 2*(diag(Cani) - Cani);
                end
            end
        end
    end
end

    function out = mcDonald(x)
        out = x.^(5/6.).*besselk(5./6,x)./(2^(5/6.)*gamma(11/6.)) ;
        out(find(x == 0)) = 3/5.;
        


