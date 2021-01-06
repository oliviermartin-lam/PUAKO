function [xStars,yStars,fStars] = detectPeaks(img,varargin)
inputs = inputParser;
inputs.addRequired('img',@isnumeric);
inputs.addParameter('threshold',3,@isnumeric);
inputs.addParameter('mask',true(size(img)),@islogical);
inputs.addParameter('border',30,@isnumeric);
inputs.addParameter('nIter',10,@isnumeric);
inputs.addParameter('nSeg1D',2^1,@isnumeric);
inputs.addParameter('nBox',30,@isnumeric);
inputs.addParameter('pupil',[],@isnumeric);
inputs.addParameter('samp',2,@isnumeric);
inputs.addParameter('nStarsMax',1e1,@isnumeric);
inputs.addParameter('iMax',[],@isnumeric);
inputs.addParameter('flagBrightestStar',false,@islogical);
inputs.parse(img,varargin{:});

threshold   = inputs.Results.threshold;
mask        = inputs.Results.mask;
border      = inputs.Results.border;
nIter       = inputs.Results.nIter;
nSeg1D      = inputs.Results.nSeg1D;
nPixPerSeg  = size(img)/nSeg1D;
nBox        = inputs.Results.nBox;
pupil       = inputs.Results.pupil;
samp        = inputs.Results.samp;
nStarsMax   = inputs.Results.nStarsMax;
iMax        = inputs.Results.iMax;
flagBrightestStar = inputs.Results.flagBrightestStar;

% Get image statistics
im_masked = img.*mask;
[~,bg,sigma_noise] = puakoTools.getFlux(im_masked); 

% Thresholding
im_thresholded = im_masked.*(im_masked > threshold*sigma_noise + bg);

% Crop borders + Median filtering
img2 = im_thresholded(border+1:end-border,border+1:end-border);
for i=1:nIter
    img2(:) = medfilt1(img2(:));
end


if flagBrightestStar
    
    % Get the max position on filtered image
    [yStars,xStars] = find(img2 == max(img2(:)));
    xStars = xStars + border;
    yStars = yStars + border;
    
    % Get the max position on the real image
    nLoc = length(xStars);
    pVal = zeros(1,nLoc);
    for i=1:nLoc
        pVal(i) = img(yStars(i),xStars(i));
    end
    xStars = {xStars(find(pVal == max(pVal)))};
    yStars = {yStars(find(pVal == max(pVal)))};
    fStars = {img(yStars{1},xStars{1})};
else
    % segment the image - squared image
    [nY,nX]= size(img2);
    xStars = cell(1,nStarsMax);
    yStars = cell(1,nStarsMax);
    fStars = cell(1,nStarsMax);
    nseg   = 0;
    k      = 0;
    
    for kx=1:nSeg1D
        % x-axis
        idx = max(1,(kx-1)*nPixPerSeg(2)+1) : min(kx*nPixPerSeg(2),nX);
        % segment x-center
        xc_seg = (kx-1)*nPixPerSeg(2) + length(idx)/2+1;
        for ky=1:nSeg1D
            %y-axis
            idy = max(1,(ky-1)*nPixPerSeg(1)+1) : min(ky*nPixPerSeg(1),nY);
            % loop on segment
            nseg = nseg + 1;
            % segment y-center
            yc_seg = (ky-1)*nPixPerSeg(1) + length(idy)/2 + 1;
            
            % take a segment
            im_seg = img2(idy,idx);
            
            if max(im_seg(:)) > threshold*sigma_noise + bg
                % Perform detection
                k=k+1;
                
                %1\ ---- Localize the PSF peak position on the filtered image
                [y0,x0] = find(im_seg == max(im_seg(:)));
                
                %2\ ---- filter hot pixels
                idy_seg     = max(1,y0(1) - nBox/2-1) : min(length(idy),y0(1) + nBox/2);
                idx_seg     = max(1,x0(1) - nBox/2-1) : min(length(idx),x0(1) + nBox/2);
                im_box      = im_seg(idy_seg,idx_seg);
                
                for jj=1:nIter
                    for n = 1:length(y0)
                        if im_seg(y0(n),x0(n)) > mean(im_box(:))*threshold
                            im_seg(y0(n),x0(n)) = 0;
                        end
                    end
                    [y0,x0] = find(im_seg == max(im_seg(:)));
                end
                
                
                if max(im_seg(:)) > threshold*sigma_noise
                    % restart from clean image
                    [y0,x0]  = find(im_seg == max(im_seg(:)));
                    
                    %3\ ---- Localize the PSF peak on the real image
                    nLoc = length(x0);
                    pVal = zeros(1,nLoc);
                    for i=1:nLoc
                        pVal(i) = img(y0(i),x0(i));
                    end
                    
                    %4\ ---- Measure the FWHM
                    [Fx,Fy,~,im_fit] = puakoTools.getFWHM(im_box,'pixelScale',1,'rebin',1,'method','Moffat','iMax',iMax);
                    
                    if abs(Fx-Fy)/Fx < 0.2  || max(im_box(:)) > iMax % single star
                        xStars{k} = x0(find(pVal == max(pVal))) - length(idx)/2  - 1 + xc_seg + border;
                        yStars{k} = y0(find(pVal == max(pVal))) - length(idy)/2  - 1 + yc_seg + border;
                        fStars{k} = img(yStars{k},xStars{k});
                    else  % binary stars
                        nS = 2;
                        [Fx,Fy,~,~,beta] = puakoTools.getFWHM(im_box,'pixelScale',1,'rebin',1,'method','Moffat','nStars',nS,'flagSymetric',true,'iMax',iMax);
                        xStars{k} = round([x0(find(pVal == max(pVal))) - length(idx)/2 - 1 + xc_seg + border + beta(nS+3:2*nS+2)]);
                        yStars{k} = round([y0(find(pVal == max(pVal))) - length(idy)/2 - 1  + yc_seg + border + beta(2*nS+3:3*nS+2)]);
                        fStars{k} = diag(img(yStars{k},xStars{k}));
                    end
                end
            end
        end
    end
    
    % concatenate results
    idGood = find(~cellfun(@isempty,xStars));
    xStars  = xStars(idGood);
    yStars  = yStars(idGood);
    fStars  = fStars(idGood);
end
