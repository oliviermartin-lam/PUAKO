function trs = restoreKeckTelemetry(trs)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));

%% 1\ Restore telemetry data
trsData = restore_idl('filename',trs.path_trs);

%% 2\ Get fits header and restore observing conditions
list  = trs.fitsHdr.field;
val   = trs.fitsHdr.value;
trs.aoMode = 'NGS';
if isfield(trsData,'B')
    trs.aoMode = 'LGS';
    trs.lgs.height = str2num(cell2mat(val(strcmp(list,'AOFCSALT'))));
end
trs.tel.zenith_angle = 90 - cell2mat(val(strcmp(list,'EL')));
trs.tel.airmass = 1/cos(trs.tel.zenith_angle*pi/180);

%% 3\ Get AO control loop data                 

%3.1. Get slopes in pixels unit
trs.wfs.slopes   = double(trsData.A.OFFSETCENTROID);
trs.wfs.nSl      = size(trs.wfs.slopes,1);
trs.wfs.nExp   = size(trs.wfs.slopes,2);

%3.2. Get DMs commands in OPD units
trs.dm.com = double(trsData.A.DMCOMMAND)*trs.dm.volt2meter;
trs.dm.nCom       = size(trs.dm.com,1);

%3.3. Get tip-tilt measurements and conversion into OPD
if ~isfield(trsData,'B')
    trs.tipTilt.slopes  = (double(trsData.A.RESIDUALWAVEFRONT(trs.dm.nCom+1:trs.dm.nCom+2,:))); %angle in arcsec
    trs.tipTilt.com = double(trsData.A.TTCOMMANDS);
else
    trs.tipTilt.slopes = (double(trsData.B.DTTCENTROIDS));
    trs.tipTilt.com = double(trsData.B.DTTCOMMANDS);
end
trs.tipTilt.slopes= trs.tipTilt.tilt2meter*trs.tipTilt.slopes;
trs.tipTilt.slopes= bsxfun(@minus,trs.tipTilt.slopes,mean(trs.tipTilt.slopes,2));
trs.tipTilt.com= trs.tipTilt.tilt2meter*trs.tipTilt.com;
trs.tipTilt.com = bsxfun(@minus,trs.tipTilt.com,mean(trs.tipTilt.com,2));
trs.tipTilt.nExp = size(trs.tipTilt.slopes,2);

%% 4\ Get system matrices and reconstructed wavefront

%4.1\ Get DM commands reconstructors from slopes
MC              = reshape(double(trsData.RX),trs.wfs.nSl ,trs.dm.nCom+3,trsData.NREC); %command matrix
trs.mat.R    = trs.dm.volt2meter*permute(MC(:,1:trs.dm.nCom,:),[2,1,3]);
trs.mat.Rtt = trs.dm.volt2meter*permute(MC(:,trs.dm.nCom+1:trs.dm.nCom+2,:),[2,1,3]);

%4.2\ Get the reconstructed wavefront in OPD and in the actuators space
trs.rec.res    = trs.dm.volt2meter*double(trsData.A.RESIDUALWAVEFRONT(1:trs.dm.nCom,:));
trs.rec.res    = bsxfun(@minus,trs.rec.res,mean(trs.rec.res,2));
trs.rec.focus = double(trsData.A.RESIDUALWAVEFRONT(end,:));
trs.rec.focus = bsxfun(@minus,trs.rec.focus,mean(trs.rec.focus,2));

% fill vector to get 21x21 actuators
u = zeros(trs.dm.nActuators^2,trs.wfs.nExp);
u(trs.dm.validActuators,:) = trs.rec.res;
trs.rec.res = u;
u = zeros(trs.dm.nActuators^2,trs.wfs.nExp);
u(trs.dm.validActuators,:) = trs.dm.com;
trs.dm.com = u;


%% 5\ Get the loop status and model transfer function
%5.1. Delays
wssmprg           = str2num(cell2mat(val(strcmp(list,'WSSMPRG'))));
[trs.holoop.lat,trs.ttloop.lat] = estimateLoopDelay(wssmprg,trs.aoMode);
%5.2. Frequency
trs.holoop.freq = 1/(100e-9*mean(diff(trsData.A.TIMESTAMP)));
trs.ttloop.freq =trs.holoop.freq;
if isfield(trsData,'B')
    trs.ttloop.freq =  1/(100e-9*mean(diff(trsData.B.TIMESTAMP)));   
end
%5.3. RTC controller HO loop
trs.holoop.gain = double(trsData.DM_SERVO(1));
trs.holoop.tf.num=  (double(trsData.DM_SERVO(1:4))');
trs.holoop.tf.den= (double(trsData.DM_SERVO(5:end))');
trs.ttloop.gain = double(trsData.DT_SERVO(1));
trs.ttloop.tf.num=  (double(trsData.DT_SERVO(1:4))');
trs.ttloop.tf.den= (double(trsData.DT_SERVO(5:end))');

end
