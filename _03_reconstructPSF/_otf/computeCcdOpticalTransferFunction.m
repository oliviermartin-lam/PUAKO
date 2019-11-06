%{
------------HEADER-----------------
Objective          ::  Compute the CCD transfer function

INPUT VARS
 psCCD         :: The detector pixel scale in arcsec/pixel
psOTF           :: the OTF pixel scale in rad^-1/pixel
nOtf              :: the OF dimension in pixels

OUTPUT VARS
otfCCD            :: Model the CCD TF

Created by      :: O. Beltramo-Martin - ONERA/LAM
Creation date  :: 11/05/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function otfCCD = computeCcdOpticalTransferFunction(psCCD,psOTF,nOtf)
inputs = inputParser;
inputs.addRequired('psCCD',@isnumeric);
inputs.addRequired('psOTF',@isnumeric);
inputs.addRequired('nOtf',@isnumeric);
inputs.parse(psCCD,psOTF,nOtf);

%1\ Creating the angular frequencies vector in rad^-1
u1D         = (-1+1/2/nOtf:2/nOtf:1-1/2/nOtf)*psOTF*nOtf/2;
[u2Dx,u2Dy] = meshgrid(u1D);
msk         = hypot(u2Dx,u2Dy)<=1;
%2\ Computing the CCD transfer function: sinc(x) = sin(pi*x)/(pi*x)
otfCCD    = sinc(u2Dx*psCCD*constants.arcsec2radian).*sinc(u2Dy*psCCD*constants.arcsec2radian);




                