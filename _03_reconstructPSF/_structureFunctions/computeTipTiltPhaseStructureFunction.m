%{
------------HEADER-----------------
Objective          ::  Compute bi-dimensional phase structure function of the aliasing phase

INPUT VARS
 psfr          :: The psfReconstruction top-level class

OUTPUT VARS
 sf_2D             :: bi-dimensional phase structure function map (Toeplitz) of the residual tip-tilt (psfr.nOtf x psfr.nOtf)
Ctt                  :: Covariance matrix of the tip-tilt error (2 x 2)
Created by         :: O. Beltramo-Martin - ONERA/LAM
Creation date      :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function obj = computeTipTiltPhaseStructureFunction(obj)
inputs = inputParser;
inputs.addRequired('obj',@(x) isa(x,'psfReconstruction'));
inputs.parse(obj);

%1\ Discretizing the pupil domain
res  = getGridCoordinates(obj.otf.nOtf,obj.otf.nOtf,0.5);
Y = res.x2D;
X = res.y2D;
%2\ Get the tip-tilt covariance matrix
obj.cov.Ctt    = obj.trs.tipTilt.com*obj.trs.tipTilt.com'/obj.trs.tipTilt.nExp;
%3\ Get the tip-tilt covariance
Du_tt   = (2*pi/obj.trs.cam.wavelength)^2*(obj.cov.Ctt - obj.trs.res.noise(1).Cn_tt);
%4\ Rotate the X/Y axis if necessary
th =  -obj.trs.tel.pupilAngle*pi/180;
Xr = X*cos(th) +Y*sin(th);
Yr = -X*sin(th) + Y*cos(th);
%5\Get the tip-tilt sf
obj.sf.Dtt   = Du_tt(1,1)*Xr.^2 + Du_tt(2,2)*Yr.^2 + Du_tt(1,2)*Xr.*Yr +Du_tt(2,1).*Yr'.*Xr';