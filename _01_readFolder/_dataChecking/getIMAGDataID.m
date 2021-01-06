function flagStatus = getIMAGDataID(obj,varargin)
inputs = inputParser;
inputs.addRequired('obj', @(x) isa(x,'puako'));
inputs.addParameter('path_imag', obj.path_imag,@ischar);
inputs.addParameter('path_dark', obj.data_folder_path.dark,@ischar);
inputs.addParameter('path_sky', obj.data_folder_path.sky,@ischar);
inputs.addParameter('path_div', obj.data_folder_path.div,@ischar);
inputs.parse(obj,varargin{:});
flagStatus = 'ok';

%% 1 Get ID of fits file data within the IMAG folder

% Check the folder
path_img = inputs.Results.path_imag;
path_dk = inputs.Results.path_dark;
path_sk = inputs.Results.path_sky;
path_d = inputs.Results.path_div;

% Read files inside
imag_list = dir(path_img);

% Get .fits file only
nDark = 0;
nSky = 0;
nDiv = 0;

fits_list = imag_list(contains({imag_list.name},'.fits'));

%% DARK FILES

if isempty(fits_list)
    idx = contains(upper({imag_list.name}),'DARK');
    if nnz(idx)
        path_dark = [path_img,imag_list(find(idx)).name,'/'];
        folder_dark = dir(path_dark);
        fits_list = folder_dark(contains({folder_dark.name},'.fits'));
        nDark = numel(fits_list);
    elseif ~isempty(path_dk)
        folder_dark = dir(path_dk);
        fits_list = folder_dark(contains({folder_dark.name},'.fits'));
        nDark = numel(fits_list);
    end
    
    %% SKY FILES
    idx = contains(upper({imag_list.name}),'SKY');
    if nnz(idx)
        path_sky = [path_img,imag_list(find(idx)).name,'/'];
        folder_sky = dir(path_sky);
        fits_list = folder_sky(contains({folder_sky.name},'.fits'));
        nSky = numel(fits_list);
    elseif ~isempty(path_sk)
        folder_sky = dir(path_sk);
        fits_list = folder_sky(contains({folder_sky.name},'.fits'));
        nSky = numel(fits_list);
    end
    
    %% DIVERSITY FILES
    idx = contains(upper({imag_list.name}),'DIV');
    if nnz(idx)
        path_div = [path_img,imag_list(find(idx)).name,'/'];
        folder_div = dir(path_div);
        fits_list = folder_div(contains({folder_div.name},'.fits'));
        nDiv = numel(fits_list);
    elseif ~isempty(path_d)
        folder_div = dir(path_d);
        fits_list = folder_div(contains({folder_div.name},'.fits'));
        nDiv = numel(fits_list);
    end
    
    uiwait(warndlg(['There is no .fits file in the folder, but I found ',num2str(nDark),' dark files and ', num2str(nSky), ' sky and ', num2str(nDiv),' files. Please check that']));
    flagStatus = 'warning';
    obj.data_folder_id.imag = [];
else
    [obj.data_folder_id.imag,flagStatus]  = grabDataID(path_img,fits_list);
end

%% 2 Check the Dark folder

idx = contains(upper({imag_list.name}),'DARK');
if nnz(idx)
    path_dark = [path_img,imag_list(find(idx)).name,'/'];
    folder_dark = dir(path_dark);
    fits_list = folder_dark(contains({folder_dark.name},'.fits'));
    [obj.data_folder_id.dark,flagStatus]  = grabDataID(path_dark,fits_list);
elseif ~isempty(path_dk)
    folder_dark = dir(path_dk);
    fits_list = folder_dark(contains({folder_dark.name},'.fits'));
    [obj.data_folder_id.dark,flagStatus]  = grabDataID(path_dk,fits_list);
end

%% 2 Check the Sky folder

idx = contains(upper({imag_list.name}),'SKY');
if nnz(idx)
    path_sky = [path_img,imag_list(find(idx)).name,'/'];
    folder_sky = dir(path_sky);
    fits_list = folder_sky(contains({folder_sky.name},'.fits'));
    [obj.data_folder_id.sky,flagStatus]  = grabDataID(path_sky,fits_list);
elseif ~isempty(path_sk)
    folder_sky = dir(path_sk);
    fits_list = folder_sky(contains({folder_sky.name},'.fits'));
    [obj.data_folder_id.sky,flagStatus]  = grabDataID(path_sk,fits_list);
end

%% 3 Check the Div folder

idx = contains(upper({imag_list.name}),'DIV');
if nnz(idx)
    path_div = [path_img,imag_list(find(idx)).name,'/'];
    folder_div = dir(path_div);
    fits_list = folder_div(contains({folder_div.name},'.fits'));
    [obj.data_folder_id.div,flagStatus]  = grabDataID(path_div,fits_list);
elseif ~isempty(path_d)
    folder_div = dir(path_d);
    fits_list = folder_div(contains({folder_div.name},'.fits'));
    [obj.data_folder_id.div,flagStatus]  = grabDataID(path_d,fits_list);
end

function [out,flagStatus] = grabDataID(path,list)

 if isempty(list)
     flagStatus = 'warning';
     out = [];
     return
 end

flagStatus = 'ok';
nIMAG = numel(list);

% Get data id and header
path_imag = cell(1,nIMAG);
name_imag  = cell(1,nIMAG);
type_imag = cell(1,nIMAG);
imagBytes = cell(1,nIMAG);
hdr = cell(1,nIMAG);
hdr_name = cell(1,nIMAG);

for k=1:nIMAG
    % Get the path
    path_imag{k} = [path,list(k).name];
    % Get the name and type
    tmp = split(list(k).name,'.');
    name_imag{k} = [tmp{1:end-1}];
    type_imag{k} = tmp{end};
    % Get the size in Mega Bytes
    imagBytes{k} = double(list(k).bytes)/1000^2;
    % Get the header    
    tmp = fitsinfo(path_imag{k});
    if isfield(tmp,'PrimaryData')
        tmp = tmp.PrimaryData.Keywords;        
        hdr{k} = tmp;
        % Get the name in the header    
        hdrname_tmp = split(cell2mat(tmp(strcmp(tmp(:,1),'FILENAME'),2)),'.');
        hdr_name{k} = hdrname_tmp{1};
    else
        fprintf('\n');
        warning(sprintf('Looks like the file %s format is ackward and there is no PrimaryData field, sorry !',tmp));
        flagStatus = 'warning';
    end
end

out  = struct('name',name_imag,'origin_name',hdr_name,'type',type_imag,'size',imagBytes,'hdr',hdr,'path',path_imag);