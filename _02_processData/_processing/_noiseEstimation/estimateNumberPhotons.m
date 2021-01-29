%{
------------HEADER-----------------
Objective          ::  Estimate the number of photons

INPUT VARS
 trs          :: The telemetry top-level class

OUTPUT VARS
trs.wfs.nph and trs.tipTilt.nph are calculated

Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 02/10/2020
                      
Change Record:     ::
------------HEADER END----------------
%}

function trs = estimateNumberPhotons(trs)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));
inputs.parse(trs);

% High order WFS
map = trs.wfs.intensity;

if any(map(:))
    % Flux estimation
    pixInt = sort(map,1);
    trs.wfs.nSl_c = nnz(trs.mat.R(100,:,1));
    pixCog = pixInt(end-trs.wfs.nSl_c/2+1:end,:);
    d = trs.dm.pitch;
    Te = 1./trs.holoop.freq;
    trs.wfs.nph = mean(pixCog(:))/d^2/Te; % in ADU/m2/s
    trs.tipTilt.nph = trs.wfs.nph;
    
    % Read-out noise estimation
    tmp = mean(pixInt,2);
    tmp2= tmp(tmp>0);
    idx = find(tmp2 == min(tmp2));
    idx2 = find(tmp == tmp2(idx));
    trs.wfs.ron  = std(pixInt(idx2,:)); %in ADU
    
    if strcmp(trs.aoMode,'LGS')
        trs.tipTilt.nph = mean(trs.tipTilt.intensity(:))/trs.tel.Dcircle*trs.ttloop.freq; % in ADU/m2/s
        trs.tipTilt.ron =  std(trs.tipTilt.intensity(:));
    end
end
%To get number of photons, I miss the detector gain and quantum
%efficiencyas well as the optical path throughput. Or a zeropoint
%calibration.