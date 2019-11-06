%{
------------HEADER-----------------
Objective          ::  Display the Ao transfer function
reconstructed phase
INPUT VARS
trs         :: the telemetry class
jindex            :: The Noll's j-indexing of Zernike modes

OUTPUT VARS
std_n            :: The standard-deviation value per mode

Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 11/01/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function displayTransferFunction(trs)
inputs = inputParser;
inputs.addRequired('trs', @(x) isa(x,'telemetry'));

% Instantiation
nObj = numel(trs);
%color
c1 = [0 0 255]/255;%black
c2 = [255,0,0]/255; %red
colors_p = [linspace(c1(1),c2(1),nObj)',linspace(c1(2),c2(2),nObj)',linspace(c1(3),c2(3),nObj)'];

h = figure;
for kObj=1:nObj
    subplot(2,2,1)
    loglog(trs(kObj).holoop.tf.freq,abs(trs(kObj).holoop.tf.ol),'Color',colors_p(kObj,:));
    hold on;
    subplot(2,2,2)
    loglog(trs(kObj).holoop.tf.freq,abs(trs(kObj).holoop.tf.ctf),'Color',colors_p(kObj,:));
    hold on;
    subplot(2,2,3)
    loglog(trs(kObj).holoop.tf.freq,abs(trs(kObj).holoop.tf.rtf),'Color',colors_p(kObj,:));
    hold on;
    subplot(2,2,4)
    loglog(trs(kObj).holoop.tf.freq,abs(trs(kObj).holoop.tf.ntf),'Color',colors_p(kObj,:));
    hold on;
end

subplot(2,2,1)
ylabel('Magnitude','interpreter','latex','FontSize',20);
xlabel('Frequency (Hz) ','interpreter','latex','FontSize',20);
legend([trs.obj_name],'interpreter','latex','FontSize',10,'Location','northeast');
title('Open-loop transfer function','interpreter','latex','FontSize',20);
grid ON
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex');

subplot(2,2,2)
ylabel('Magnitude','interpreter','latex','FontSize',20);
xlabel('Frequency (Hz) ','interpreter','latex','FontSize',20);
legend([trs.obj_name],'interpreter','latex','FontSize',10,'Location','northwest');
title('Closed-loop transfer function','interpreter','latex','FontSize',20);
grid ON
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex');

subplot(2,2,3)
ylabel('Magnitude','interpreter','latex','FontSize',20);
xlabel('Frequency (Hz) ','interpreter','latex','FontSize',20);
legend([trs.obj_name],'interpreter','latex','FontSize',10,'Location','northwest');
title('Rejection transfer function','interpreter','latex','FontSize',20);
grid ON
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex');

subplot(2,2,4)
ylabel('Magnitude','interpreter','latex','FontSize',20);
xlabel('Frequency (Hz) ','interpreter','latex','FontSize',20);
legend([trs.obj_name],'interpreter','latex','FontSize',10,'Location','northwest');
title('Noise transfer function','interpreter','latex','FontSize',20);
grid ON
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex');