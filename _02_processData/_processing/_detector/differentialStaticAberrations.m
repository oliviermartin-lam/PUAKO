function diff_field_stat = differentialStaticAberrations(src,cam,varargin)
inputs = inputParser;
inputs.addRequired('src',@isstruct);
inputs.addRequired('cam',@isstruct);
inputs.addParameter('interpMethod','nearest',@ischar);
inputs.addParameter('nWeight',4,@isnumeric);
inputs.parse(src,cam,varargin{:});
interpMethod = inputs.Results.interpMethod;
nWeight = inputs.Results.nWeight;

%1\ Get static maps calibration
map     = cam.field_static_map;
xmap    = cam.x_stat; % in arcsec
ymap    = cam.y_stat;

%2\ Get source position in the CCD during observation
xsrc    = src(1).x(1); % in arcsec
ysrc    = src(1).y(1);

%3\ Get source position in the CCD during NCPA calibration
%xncpa = cam.x_ncpa;
%yncpa = cam.y_ncpa;

xncpa = xsrc;% so far I don't know the position
yncpa = ysrc;

%4\ Select the four closest neighbors around the src position
d_src                           = hypot(xmap - xsrc,ymap - ysrc);
[weight_src,idx_closest_src]    = sort(d_src);
map_src                         = map(:,:,idx_closest_src(1:nWeight));
weight_src                      = 1./weight_src(1:nWeight);
weight_src                      = weight_src/sum(weight_src);

%4\ Select the four closest neighbors around the stat position
d_ncpa                          = hypot(xmap - xncpa,ymap - yncpa);
[weight_ncpa,idx_closest_ncpa]  = sort(d_ncpa);
map_ncpa                        = map(:,:,idx_closest_ncpa(1:nWeight));
weight_ncpa                     = 1./weight_ncpa(1:nWeight);
weight_ncpa                     = weight_ncpa/sum(weight_ncpa);

%5\ Interpolation
if strcmp(interpMethod,'nearest')
    % weighting 
    mapSrc  = sum(reshape(map_src,[],nWeight).*weight_src,2);
    mapNcpa = sum(reshape(map_ncpa,[],nWeight).*weight_ncpa,2);
    % difference
    nRes            = sqrt(size(mapNcpa,1));
    diff_field_stat = reshape(mapSrc - mapNcpa,nRes,nRes);
else
    % No alternative so far
    diff_field_stat = 0;     
end
