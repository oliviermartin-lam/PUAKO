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

function obj = computeAnisoplanatismPhaseStructureFunction(obj)
inputs = inputParser;
inputs.addRequired('obj',@(x) isa(x,'psfReconstruction'));
inputs.parse(obj);


obj.sf.Dani_l   = 0;
obj.sf.DaniTT_l = 0;
obj.sf.Dani     = 0;
obj.otf.Kani    = 1;

%1\ Verify if there is anisoplanatism
xSrc                        = cell2mat({obj.trs.src.x});
ySrc                        = cell2mat({obj.trs.src.y});
obj.flags.isFocalAniso      = strcmp(obj.trs.aoMode,'LGS');
if obj.flags.isFocalAniso
    obj.flags.isAngularAniso = any(~isequal(obj.trs.lgs.x,xSrc)) || any(~isequal(obj.trs.lgs.y,ySrc));
    obj.flags.isTiltAniso    = any(~isequal(obj.trs.ngs.x,xSrc)) || any(~isequal(obj.trs.ngs.y,ySrc));
else
    obj.flags.isAngularAniso = any(~isequal(obj.trs.ngs.x,xSrc)) || any(~isequal(obj.trs.ngs.y,ySrc));
    obj.flags.isTiltAniso    = false;
end
obj.flags.isAniso           = obj.flags.isFocalAniso + obj.flags.isAngularAniso + obj.flags.isTiltAniso;

if obj.flags.isAniso    
    r053 =  obj.trs.res.seeing.r0^(-5/3)  *(obj.trs.atm.wavelength/obj.trs.cam.wavelength)^2;
    %Note: I do consider that the seeing estimated from the telemetry applies to all layers
    Cn2 = obj.trs.atm.weights*r053;
    
    if obj.flags.isAngularAniso && ~obj.flags.isFocalAniso
        %1\ Get the angular-anisoplanatism phase structure function
        obj.sf.Dani_l = instantiateAnisoplanatism(obj,obj.trs.ngs);   
    
    elseif obj.flags.isFocalAniso
        %2\ Get the focal-angular-anisoplanatism phase structure function
        obj.sf.Dani_l     = instantiateAnisoplanatism(obj,obj.trs.lgs);
        lgsInf            = obj.trs.lgs;
        % get the angular anisoplanatism contribution
        lgsInf.height     = Inf;
        tmp               = instantiateAnisoplanatism(obj,lgsInf);
        obj.otf.KaniAng   = exp(-0.5*squeeze(sum(bsxfun(@times, tmp , reshape(Cn2,1,1,[])), 3)));
    end
    
    obj.sf.Dani = squeeze(sum(bsxfun(@times, obj.sf.Dani_l , reshape(Cn2,1,1,[])), 3));
    
    if obj.flags.isTiltAniso
        %3\ Get the tip-tilt anisoplanatism phase structure function
        obj.sf.DaniTT_l     = instantiateAnisoplanatism(obj,obj.trs.ngs,'isTT',true);
        obj.sf.DaniTT       = sum(bsxfun(@times, obj.sf.DaniTT_l, reshape(Cn2,1,1,[])), 3);
        obj.otf.KaniTT      = exp(-0.5*obj.sf.DaniTT);
        obj.sf.Dani         = squeeze(obj.sf.Dani+ obj.sf.DaniTT);
        obj.sf.Dani_l       = obj.sf.Dani_l  + obj.sf.DaniTT_l;
    end
    
    %4\ Get the anisoplanatism spatial filter
    if obj.flags.toeplitz
        obj.otf.Kani                = exp(-0.5*obj.sf.Dani);
    else
        idx                         = ones(obj.otf.dk);
        dp                          = obj.trs.tel.Ddm/(obj.otf.dk-1);
        obj.otf.otfDen              = modes2Otf(obj.sf.Dani*0,obj.trs.mat.dmIF',idx,obj.otf.nOtf,'D',obj.trs.tel.Ddm,'dk',dp,'method','zonal');
        obj.otf.mskOtf              = obj.otf.otfDen > 1e-6;
        Kani                        = modes2Otf(obj.sf.Dani,obj.trs.mat.dmIF',idx,obj.otf.nOtf,'D',obj.trs.tel.Ddm,'dk',dp,'method','zonal');
        obj.otf.Kani                = 0*Kani;
        obj.otf.Kani(obj.otf.mskOtf)   = Kani(obj.otf.mskOtf)./obj.otf.otfDen(obj.otf.mskOtf);
    end
end
