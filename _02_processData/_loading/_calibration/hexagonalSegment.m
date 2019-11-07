function out = hexagonalSegment(X0,Y0,segsize,rotAngle,varargin)
inputs = inputParser;
inputs.addRequired('X0', @isnumeric);
inputs.addRequired('Y0', @isnumeric);
inputs.addRequired('segsize', @isnumeric);
inputs.addRequired('rotAngle', @isnumeric);
inputs.addParameter('overSamp', 2,@isnumeric);
inputs.parse(X0,Y0,segsize,rotAngle,varargin{:});

overSamp = inputs.Results.overSamp;
%1\ Include rotation
th = rotAngle*pi/180;
    
if  any([0,90,180,270,360,-90,-180,-270,-360] == th)
    out = (abs(X0) <= sqrt(3)*segsize/4).*(abs(Y0) <= X0/sqrt(3) + segsize/2) .*(abs(Y0) <= -X0/sqrt(3) + segsize/2);
    out = tools.rotateIm(out,th);
else
    %1 interpolate the dimension to mitigate pixellization effects
    overSamp = max(round(400/size(X0,1)),1); %% 400 pixels is enough to have a good description
    res = getGridCoordinates(size(X0,1)*overSamp,size(Y0,1)*overSamp,max(abs(X0(:))));    
    %2\ Hexagonal shape with rotation
    Xr = res.x2D.*cos(th) + res.y2D.*sin(th);
    Yr = res.y2D.*cos(th) - res.x2D.*sin(th);
    out = (abs(Xr) <= sqrt(3)*segsize/4).*(abs(Yr) <= Xr/sqrt(3) + segsize/2) .*(abs(Yr) <= -Xr/sqrt(3) + segsize/2);
    %3\ Get back to the original resolution
    out = tools.interpolate(out,[size(X0,1),size(Y0,1)]);
end