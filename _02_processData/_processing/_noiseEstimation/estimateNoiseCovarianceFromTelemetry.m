%{
------------HEADER-----------------
Objective          ::  Estimate the noise covariance matrix

INPUT VARS
 trs          :: The telemetry top-level class
method    :: (optional) 'method': 'autocorrelation', 'interpolation' or 'rtf'

OUTPUT VARS
Cn                   :: The HO WFS noise covariance matrix (psfr.nActu^2 x psfr.nActu^2)
Cn_tt               :: The TT WFS noise covariance matrix (2 x 2)
var_n               :: The HO WFS noise variance in meter^2
var_n_tt           :: The TT WFS noise variance in meter^2
Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function res = estimateNoiseCovarianceFromTelemetry(trs,varargin)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));
inputs.addParameter('flagNoisemethod','autocorrelation',@ischar);
inputs.parse(trs,varargin{:});
method = inputs.Results.flagNoisemethod;

res = [];

if strcmp(method,'theoretical')
    % !!!! to be review !!!1
    %Number of photons
    
    d       = trs.dm.pitch;
    ron     = trs.wfs.ron;
    wvl     = trs.wfs.wavelength;
    psInRad = constants.arcsec2radian*trs.wfs.pixelScale/1e3;
    
    % Estimation of the number of photon
    nph = trs.wfs.nph*trs.dm.pitch^2/trs.holoop.freq;
    if strcmp(trs.aoMode,'NGS')
        nph_tt = nph;
        ron_tt = ron;
        d_tt = d;
        psInRad_tt = psInRad;
        wvl_tt = wvl;
    end
    
    % Noise variance in pix^2 between
    varPix = (4*pi^2*(ron/nph)^2 + pi^2*1./nph )*(wvl/(2*pi*d*psInRad))^2; %nph in ADU... not ok I think
    %Propagation through the reconstructor    
    res.Cn_ho = varPix*trs.mat.R*trs.mat.R';    
    validInput = find(trs.mat.R(:,100)~=0);
    
    varPix_tt = (4*pi^2*(ron_tt/nph_tt)^2 + pi^2*1./nph_tt )*(wvl_tt/(2*pi*d_tt*psInRad_tt))^2;
    res.Cn_tt = varPix_tt*trs.mat.Rtt*trs.mat.Rtt';

elseif strcmp(method,'nonoise')
    res.Cn_ho = zeros(size(trs.dm.com,1));
    res.Cn_tt = zeros(size(trs.tipTilt.com,1));
    validInput = std(trs.dm.com,[],2)~=0;

else 
    %1\ HO WFS noise covariance matrix
    res.Cn_ho = getNoiseCovariance(trs.dm.com,'flagNoisemethod',method,'rtf',trs.holoop.tf.ctf./trs.holoop.tf.wfs,'FSAMP',trs.holoop.freq);
    validInput = std(trs.dm.com,[],2)~=0;

    %2\ TT WFS noise covariance matrix
    res.Cn_tt = getNoiseCovariance(trs.tipTilt.com,'flagNoisemethod',method,'rtf',trs.ttloop.tf.ctf./trs.ttloop.tf.wfs,'FSAMP',trs.ttloop.freq);
end

res.std_ho = sqrt(trace(res.Cn_ho(validInput,validInput))/nnz(validInput))*1e9;
res.std_tt = sqrt(trace(res.Cn_tt))*1e9;
res.method = method;


function Cnn = getNoiseCovariance(s,varargin)
inputs = inputParser;
inputs.addRequired('s',@isnumeric);
inputs.addParameter('FSAMP',1e3,@isnumeric);
inputs.addParameter('flagNoisemethod','autocorrelation',@ischar);
inputs.addParameter('nfit',1,@isnumeric);
inputs.addParameter('nshift',1,@isnumeric);
inputs.addParameter('rtf',1,@isnumeric);
inputs.parse(s,varargin{:});

method = inputs.Results.flagNoisemethod;
nfit = inputs.Results.nfit;
nshift = inputs.Results.nshift;
rtf = inputs.Results.rtf;

%1\ Mean removal
s       = squeeze(s);
s       = bsxfun(@minus,s,mean(s,2));
validInput = std(s,[],2)~=0;
[nS,nF] = size(s);
Cnn     = zeros(nS);

%2\ Covariance derivation
if strcmp(method,'interpolation') % Fitting the first samples of the auto-correlation function with a polynomial model
    delay   = linspace(0,1,nfit+2);
    for i=1:nS
        g      = (ifft(fft(s(i,:)).*conj(fft(s(i,:))))/nF);
        mx     = max(g(:));
        fit    = polyfit(delay(2:end),g(2:nfit+2),nfit);
        yfit   = 0*delay;
        for k =1:nfit+1
            yfit   = yfit + fit(k)*delay.^(nfit-k+1);
        end
        Cnn(i,i) = mx - yfit(1);
    end
    Cnn = diag(Cnn);
    
elseif strcmp(method,'autocorrelation_pos') % Deriving the temporal cross-correlation difference over 2xnshift
    % If the turbulence is frozen over the  WFS frame rate, the difference of the autocorrelation and the 1-frame shifted correlation
    % is the noise variance.    
    ds_p  = s - circshift(s,[0,nshift]);
    Cnn = s*ds_p'/nF;    
    
elseif strcmp(method,'autocorrelation') % Deriving the temporal cross-correlation difference over 2xnshift
    % If the turbulence is frozen over the  WFS frame rate, the difference of the autocorrelation and the 1-frame shifted correlation
    % is the noise variance.
    ds_n  = s - circshift(s,[0,-nshift]);
    ds_p  = s - circshift(s,[0,nshift]);
    Cnn = 0.5*(s*ds_n' + s*ds_p')/nF;
    
elseif strcmp(method,'rtf') % Adjusting the noise level through the noise rejection transfer function model
    fftS = fft(s,[],2); % we verified that std(s(k,:))^2 = sum(abs(fftS(k,:)).^2): Parseval
    fftS = fftS(:,1:floor(end/2))./rtf;
    % truncate to keep high-spatial frequencies only
    sub_dim = floor(9*nF/2/10);
    fftS = fftS(:,sub_dim:end);
    cfftS= conj(fftS);
       
    for i=1:nS 
        if validInput(i)
            % Get the cross PSD
            crossPSD  = bsxfun(@times,fftS(i,:),cfftS(i:end,:));
            % Estimate the noise plateau
            Cnn(i,i:end)  = median(real(crossPSD),2);
        end
    end
    Cnn = transpose(Cnn) + Cnn - diag(diag(Cnn));    
    % Normalization
    Cnn = Cnn/nF;
end
