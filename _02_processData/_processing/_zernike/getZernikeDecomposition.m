%{
------------HEADER-----------------
Objective          ::  Get the Zernike coefficients rms value of the
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

function std_n = getZernikeDecomposition(trs,jindex)
inputs = inputParser;
inputs.addRequired('trs', @(x) isa(x,'telemetry'));
inputs.addRequired('jindex', @isnumeric);

% Instantiation
nObj = numel(trs);
jindex = jindex(jindex~=1);
nJ = length(jindex);
std_n = zeros(nObj,nJ);
%color
c1 = [0 0 255]/255;%black
c2 = [255,0,0]/255; %red
%c2 = [255,255,0]/255;% yellow
%c2 = [255,128,0]/255;% orange
colors_p = [linspace(c1(1),c2(1),nObj)',linspace(c1(2),c2(2),nObj)',linspace(c1(3),c2(3),nObj)'];

% Define the Zernike modes
nPup = sqrt(size(trs(1).mat.dmIF,1));
z = zernike(jindex(jindex~=2 & jindex~=3),nPup);
u2z      = pinv(full(z.modes))*trs(1).mat.dmIF;

for kObj=1:nObj
    % Projection onto Zernike
    z_coefs = u2z*trs(kObj).rec.res;
    tmp = std(1e9*z_coefs,[],2)';
    if find(jindex==2) || find(jindex==3)
        std_tt = 1e9*std(trs(kObj).tipTilt.com,[],2)';
        std_n(kObj,:) = [std_tt tmp];
    else
        std_n(kObj,:) = tmp;
    end
end

% Display
h = figure;
for kObj=1:nObj
    semilogy(jindex,std_n(kObj,:),'LineStyle','--','Marker','s','Color',colors_p(kObj,:),'MarkerFaceColor',colors_p(kObj,:),'MarkerSize',5);
    hold on;
end
ylabel('Zernike rms value (nm)','interpreter','latex','FontSize',20);
xlabel('j-index ','interpreter','latex','FontSize',20);
legend([{trs.obj_name}],'interpreter','latex','FontSize',10,'Location','northeast');
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex');
