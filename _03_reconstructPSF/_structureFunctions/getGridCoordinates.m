%{
------------HEADER-----------------
Objective          ::  Define coordinates in any domain.

INPUT VARS
 nX,nY          :: Dimension in pixels
pScale          :: pixel scale

OUTPUT VARS
res             :: Structure containing the cartesian (res.x,res.y) and
polar(res.rho,res.theta) coordinates

Created by         :: O. Beltramo-Martin - ONERA/LAM
Creation date      :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function [x2D,y2D,r2D,th2D,x1D,y1D] = getGridCoordinates(nX,nY,pScale,varargin)
inputs = inputParser;
inputs.addRequired('nX', @isnumeric);
inputs.addRequired('nY', @isnumeric);
inputs.addRequired('pScale', @isnumeric);
inputs.parse(nX,nY,pScale,varargin{:});


x1D = linspace(-1,1-iseven(nX)*2/nX,nX)*pScale;
y1D = linspace(-1,1-iseven(nY)*2/nY,nY)*pScale;        
[x2D,y2D] = meshgrid(x1D,y1D);
[th2D,r2D] = cart2pol(x2D,y2D);
