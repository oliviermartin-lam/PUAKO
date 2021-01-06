%{
------------HEADER-----------------
Objective          ::  Select the valid actuators regarding the outer and
inner pupils size
INPUT VARS
uin              :: the actuators commands
D1,D2          :: Outer and inner pupil diameter in meter
nX,nY          :: (optional) Number of 1D actuators over the full pupil x/y-axis
pitch          :: (optional)  DM actuators pitch

OUTPUT VARS
uout            :: The truncated actuators command vector
Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 11/01/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function [uout,mskModes] = selectActuators(uin,D1,D2,varargin)
inputs = inputParser;
inputs.addRequired('wvf',@isnumeric);
inputs.addRequired('D1',@isnumeric);
inputs.addRequired('D2',@isnumeric);
inputs.addParameter('nX',21,@isnumeric);
inputs.addParameter('nY',21,@isnumeric);
inputs.addParameter('pitch',0.5625,@isnumeric);
inputs.parse(uin,D1,D2,varargin{:});

%1\ Define the actuators grid
nX = inputs.Results.nX;
nY = inputs.Results.nY;
actuPitch = inputs.Results.pitch;
pScale = ((nX-1)+0.54)*actuPitch/2;
%1D DM actuators position for which x1D(nX/2+1) = 0 and y1D(nY/2+1) = 0;
res = getGridCoordinates(nX,nY,pScale);
r2D = res.r2D;

%2\ Define DM actuators masks
mskModes = r2D < D1/2 & r2D > D2/2;
uout = bsxfun(@times,uin,mskModes(:));