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

function res = getGridCoordinates(nX,nY,pScale,varargin)
inputs = inputParser;
inputs.addRequired('nX', @isnumeric);
inputs.addRequired('nY', @isnumeric);
inputs.addRequired('pScale', @isnumeric);
inputs.parse(nX,nY,pScale,varargin{:});

%1D
res.x1D = linspace(-1,1-~mod(nX,2)*2/nX,nX)*pScale;
res.y1D = linspace(-1,1-~mod(nY,2)*2/nY,nY)*pScale;        
%2D
[res.x2D,res.y2D] = meshgrid(res.x1D,res.y1D);

if nX == nY
    [res.th1D,res.r1D] = cart2pol(res.x1D ,res.y1D);
    [res.th2D,res.r2D] = cart2pol(res.x2D,res.y2D);
end
