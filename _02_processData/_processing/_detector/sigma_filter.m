function out = sigma_filter(image,varargin)
inputs = inputParser;
inputs.addRequired('image',@isnumeric);
inputs.addParameter('box_width',3,@isnumeric);
inputs.addParameter('nSigma',3,@isnumeric);
inputs.addParameter('flagAllPixels',false,@islogical);
inputs.addParameter('flagIterate',true,@islogical);
inputs.parse(image,varargin{:});
flagIterate = inputs.Results.flagIterate;
nSigma = inputs.Results.nSigma;

out = image;
if flagIterate && (nSigma > 2)
    for k=1:10
        out = subProcess(out,varargin{:});
    end
else
    out = subProcess(out,varargin{:});
end


function out = subProcess(image,varargin)
inputs = inputParser;
inputs.addRequired('image',@isnumeric);
inputs.addParameter('box_width',3,@isnumeric);
inputs.addParameter('nSigma',3,@isnumeric);
inputs.addParameter('flagAllPixels',false,@islogical);
inputs.addParameter('flagIterate',true,@islogical);
inputs.parse(image,varargin{:});

box_width = 2*(fix(inputs.Results.box_width)/2) + 1;
flagAllPixels = inputs.Results.flagAllPixels;

if box_width < 3    
    return
end
bw2 = box_width^2;
out=( filter_image( image,'box_width',box_width,'flagAllPixels', flagAllPixels)*bw2 - image )/(bw2-1);

nSigma = inputs.Results.nSigma;
if numel(nSigma)~=1
    nSigma = 3;
end

if nSigma == 0
    return
end

imdev = (image - out).^2;
fact = floor( nSigma^2 )/(bw2-2);
imvar = fact*( filter_image( imdev,'box_width',box_width,'flagAllPixels',flagAllPixels)*bw2 - imdev );

wok = find( imdev < imvar);
nok = nnz(wok);
npix = length( image(:) );
if nok == npix
    out = image;
    return
end
if  nok > 0
    out(wok) = image(wok);
end