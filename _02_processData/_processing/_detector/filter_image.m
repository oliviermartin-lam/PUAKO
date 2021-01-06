function out = filter_image(image,varargin)
inputs = inputParser;
inputs.addRequired('image',@isnumeric);
inputs.addParameter('box_width',3,@isnumeric);
inputs.addParameter('nSigma',3,@isnumeric);
inputs.addParameter('flagAllPixels',true,@islogical);
inputs.addParameter('flagIterate',true,@islogical);
inputs.parse(image,varargin{:});

out = image;
for i=1:(inputs.Results.box_width)/2 %(width_smooth>3)/2
    out = subProcess( out, varargin{:});
end
   

function out = subProcess(image,varargin)
inputs = inputParser;
inputs.addRequired('image',@isnumeric);
inputs.addParameter('box_width',3,@isnumeric);
inputs.addParameter('nSigma',3,@isnumeric);
inputs.addParameter('flagAllPixels',false,@islogical);
inputs.addParameter('flagIterate',true,@islogical);
inputs.parse(image,varargin{:});

width_smooth = inputs.Results.box_width;
flagAllPixels = inputs.Results.flagAllPixels;

out = image;
sim = size( image );
Lx = sim(1);
Ly = sim(2);


if inputs.Results.flagIterate
    if numel( width_smooth ) ~=1
        return;
    end
    if width_smooth < 1
        return;
    end
end

    box_wid = width_smooth;% > 3;
    
    if box_wid <3
        return
    end
    
    if flagAllPixels
        box_wid = fix( box_wid );
        radius = (box_wid/2) > 1; %???
        Lxr = Lx+radius +1;
        Lyr = Ly+radius+1;
        rr = 2*radius;
        imf = zeros( sim(1)+rr, sim(2)+rr );
        imf(radius,radius) = image;  %
        imf(1,1) = rotateImage( imf(radius:rr,:), 90 ) ;     %Left
        imf(Lxr,1) = rotateImage( imf(Lx:Lxr,:), -90 )  ;       %right
        imf(1,1) = rotateImage( imf(:,radius:rr), 180 );      %bottom
        imf(1,Lyr) = rotateImage( imf(:,Ly:Lyr), -180 ) ;        %top
    else
        radius = 0;
        imf = image;
    end
    
    out(:) = smooth( imf(:), width_smooth);
    
    if radius>0
        out = imf(radius:(Lx+radius), radius:(Ly+radius));
    end

