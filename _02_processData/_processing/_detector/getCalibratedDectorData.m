function out = getCalibratedDectorData(path_calib,headerImage,headerFields,varargin)
inputs = inputParser;
inputs.addRequired('path_calib',@ischar);          
inputs.addRequired('headerImage',@iscell); 
inputs.addRequired('headerFields',@iscell); 
inputs.addParameter('average',true,@islogical); 
inputs.parse(path_calib,headerImage,headerFields,varargin{:});

fAvg = inputs.Results.average;
% List the folders and files name here
folder_calib = dir(path_calib);

% read all the fits header
idx = find(contains(upper({folder_calib.name}),'.FITS'));
calib_name = {folder_calib(idx).name};
nFits = numel(calib_name);

% grab the calibration parametrization
nHdr = numel(headerFields);
tmp_hdr = cell(nHdr,nFits);
tmp_img = cell(nHdr,1);
idx_good = ones(1,nFits);
nPix        = zeros(2,nFits);

nPix_im = [cell2mat(headerImage(find(strcmp(headerImage(:,1),'NAXIS1')),2)),cell2mat(headerImage(find(strcmp(headerImage(:,1),'NAXIS2')),2))];
        
for j=1:nHdr
    tmp_img{j} = cell2mat(headerImage(find(strcmp(headerImage(:,1),headerFields{j})),2));
    for k=1:nFits        
        % Read keywords
        hdr = fitsinfo([path_calib,folder_calib(idx(k)).name]);
        hdr = hdr.PrimaryData.Keywords;
        tmp_hdr{j,k} = cell2mat(hdr(find(strcmp(hdr(:,1),headerFields{j})),2));        
        % Get the image size
        nPix(1,k) = cell2mat(hdr(find(strcmp(hdr(:,1),'NAXIS1')),2));        
        nPix(2,k) = cell2mat(hdr(find(strcmp(hdr(:,1),'NAXIS2')),2));        
        %check if it intersects with observing settings
        idx_good(k) = idx_good(k).*isequal(tmp_hdr{j,k},tmp_img{j}).*(nPix(1,k) == nPix_im(1)).*(nPix(2,k) == nPix_im(2));
    end
end

%Select and average the good fits files
if any(idx_good)
    id_good = idx(find(idx_good));
    nGood = nnz(idx_good);    
    if fAvg % averaging fits file
        out = zeros(nPix_im(2),nPix_im(1));
        for k=1:nGood
            out = out + fitsread([path_calib,folder_calib(id_good(k)).name])/nGood;
        end
    else % taking them all
        out = cell(1,nGood);
        for k=1:nGood
            out{k} = fitsread([path_calib,folder_calib(idx(id_good(k))).name]);
        end
    end        
else
    out = 0;
    fprintf('No valid dark calibration available for your data !\n');
end
