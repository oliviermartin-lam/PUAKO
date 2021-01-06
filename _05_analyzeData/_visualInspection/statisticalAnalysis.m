classdef statisticalAnalysis < handle
    
    properties
        camName;
        psfrFlag;
        primeFlag;
        obj_name;     
        % Folders
        flagMultipleFolders;
        paths;
        DATE;
        AOMODE;
        % SYSTEM PARAMETERS
        SYSGAIN;
        SYSFREQ;
        SYSLAT;
        SYSGAINT;
        SYSFREQT;
        SYSLATT;
        SYSAIRM;
        SYSEL;
        SYSAZ;
        SYSHLGS ;
        SYSEXP;
        SYSFOC;
        %MASS DIMM
        MASSR0;
        MASSCN2;
        MASSALT;
        % TELEMETRY PROCESSING RESULTS
        TRSR0;
        TRSDR0;
        TRSL0;
        TRSDL0;        
        TRSSEEV;
        TRSDSEEV;
        TRSSEEK;
        TRSDSEEK;
        TRSMNO;
        TRSMNOTT;
        TRSTIP;
        TRSTILT;
        TRSHO;
        TRSNPH;
        TRSNPHTT;
        % IMAGE PROCESSING RESULTS
        SKYSR;
        SKYDSR;
        SKYDFW;
        SKYFWX;
        SKYFWY;
        % IMAGE GAUSSIAN-FITTING RESULTS
        SKYGSR;
        SKYGDSR;
        SKYGFWX;
        SKYGDFWX;
        SKYGFWY;
        SKYGDFWY;
        SKYGX;
        SKYGY;
        SKYGDX;
        SKYGDY;
        SKYGF;
        SKYGDF;
        SKYGFVU;
        % IMAGE MOFFAT-FITTING RESULTS
        SKYMSR;
        SKYMDSR;
        SKYMFWX;
        SKYMDFWX;
        SKYMFWY;
        SKYMDFWY;
        SKYMX;
        SKYMY;
        SKYMDX;
        SKYMDY;
        SKYMF;
        SKYMDF;
        SKYMFVU;
        % PSF-R RESULTS
        RECSR;
        RECFWX;
        RECFWY;
        RECDFWX;
        RECMAR;
        RECPAR;
        RECWFE;
        RECNCPA;
        RECFIT;
        RECLAG;
        RECNOI;
        RECALIAS;
        RECANISO;
        RECTT;
        RECNOITT;
        RECX;
        RECY;
        RECDX;
        RECDY;
        RECF;
        RECDF;
        RECFVU;
        % PRIME RESULTS
        PRISR;
        PRIFWX;
        PRIFWY;
        PRIDFWX;
        PRIMAR;
        PRIPAR;
        PRIWFE;
        PRINCPA;
        PRIFIT;
        PRILAG;
        PRINOI;
        PRIALIAS;
        PRIANISO;
        PRITT;
        PRINOITT;
        PRIX;
        PRIY;
        PRIDX;
        PRIDY;
        PRIF;
        PRIDF;
        PRIFVU;
        PRIDGAO;
        PRIGAO;
        PRIGAOR;
        PRIDGAOR;
        PRIGTT;
        PRIDGTT;
        PRIGAL;
        PRIDGAL;
        PRIR0;
        PRIDR0;
        PRICN2;
        PRIDCN2;
        PRICOEF;
        PRIDCOEF;
    end
    
    methods
        
        function obj = statisticalAnalysis(path_save,varargin)
            inputs = inputParser;
            inputs.addRequired('path_save', @(x) isa(x,'char') || isa(x,'cell') );
            inputs.parse(path_save,varargin{:});
            
            
            % Check if folders exist            
            obj.flagMultipleFolders =  isa(path_save,'cell');
            if ~obj.flagMultipleFolders
                path_save = {path_save};
            end
            nFolders = numel(path_save);
            nFiles = 0;
            path_save_full = [];
            kk = 0;
            
            for kFolder=1:nFolders
                if ~isfolder(path_save{kFolder})
                    fprintf('Sorry, you must provide a valid folder path\n');
                    return
                else                    
                    folder_save = dir(path_save{kFolder});
                    idx = find(contains({folder_save.name},'.fits'));
                    if isempty(idx)
                        fprintf('There''s not .fits file in the folder\n');
                        return
                    else
                        kk = kk + 1;
                        path_save_full{kk} = path_save{kFolder};
                        nFiles = nFiles + numel(idx);
                    end                
                end
            end
            
            % Grab data
            nFolders = numel(path_save_full);
            obj.camName = cell(1,nFiles);
            obj.psfrFlag = zeros(1,nFiles);
            obj.primeFlag = zeros(1,nFiles);
            obj.obj_name = cell(1,nFiles);
            obj.AOMODE = cell(1,nFiles);
            obj.DATE = cell(1,nFiles);
            kk = 0;
            
            for kFolder=1:nFolders
                folder_save = dir(path_save_full{kFolder});
                idx = find(contains({folder_save.name},'.fits'));                           
                % Get .fits file name
                nG  = numel(idx);
                fName = fieldnames(obj);
                nField = length(fName)-4;                
                tmp = strsplit(path_save_full{kFolder},'/');                
                for k=1:nG
                    % Name
                    kk = kk + 1;
                    nm_tmp = split(folder_save(idx(k)).name,'_');
                    obj.DATE{kk} = tmp{end-1};
                    obj.camName{kk} = nm_tmp{1};
                    obj.obj_name{kk} = nm_tmp{2};
                    obj.primeFlag(kk) = contains(upper(nm_tmp{3}),'PRIME');
                    obj.psfrFlag(kk) = contains(upper(nm_tmp{3}),'PSFR');
                    % Fields
                    hdr = fitsinfo([path_save_full{kFolder},folder_save(idx(k)).name]);
                    hdr = hdr.PrimaryData.Keywords;
                    for j=1:nField
                        old = getfield(obj,fName{j+4});
                        new = cell2mat(hdr(strcmp(hdr(:,1),fName{j+4}),2))';
                        setfield(obj,fName{j+4},[old,new]);
                    end
                    obj.AOMODE{kk} = 'NGS';
                    if cell2mat(hdr(strcmp(hdr(:,1),'SYSHLGS'),2)) ~=0
                        obj.AOMODE{kk} = 'LGS';
                    end                    
                end
            end          
        end
        
        function displayPlot(obj,fieldx,fieldy,varargin)
            inputs = inputParser;
            inputs.addRequired('obj', @(x) isa(x,'statisticalAnalysis'));
            inputs.addRequired('fieldx', @ischar );
            inputs.addRequired('fieldy', @ischar );
            inputs.addParameter('sortby',[], @ischar );
            inputs.addParameter('resultsof','psfr', @ischar );
            inputs.addParameter('fontsize',20, @isnumeric );
            inputs.addParameter('legendx',fieldx, @ischar );
            inputs.addParameter('legendy',fieldy, @ischar );
            inputs.addParameter('xyline',false, @islogical );
            inputs.addParameter('Color','b', @ischar );
            inputs.addParameter('Marker',[], @ischar );
            inputs.addParameter('MarkerSize',7, @isnumeric );
            inputs.addParameter('thresx',[], @isnumeric );
            inputs.addParameter('thresy',[], @isnumeric );
            inputs.addParameter('scale','linear', @ischar );
            inputs.addParameter('legendLoc','northwest', @ischar );
            inputs.addParameter('polyfitOrder',0, @isnumeric );
            inputs.parse(obj,fieldx,fieldy,varargin{:});            
                                   
            %1\ Parsing inputs
            resultsof = inputs.Results.resultsof;
            fontsize = inputs.Results.fontsize;
            legendx = inputs.Results.legendx;
            legendy = inputs.Results.legendy;
            legendLoc =  inputs.Results.legendLoc;
            xyline = inputs.Results.xyline;            
            mk = inputs.Results.Marker;            
            mkc = inputs.Results.Color;            
            mks = inputs.Results.MarkerSize;            
            thresx = inputs.Results.thresx;
            thresy = inputs.Results.thresy;
            sortby = inputs.Results.sortby;
            scale = inputs.Results.scale;
            polyfitOrder = inputs.Results.polyfitOrder;
            
            if ~strcmpi(resultsof,'PSFR') && ~strcmpi(resultsof,'PRIME') && ~strcmpi(resultsof,'ALL')
                fprintf('Sorry, I can only plot PSF-R, PRIME or both results\n');
                return
            end
            
            %% Check if the fields exist
            flagerror = false;
              if ~isprop(obj,fieldx)
                  list = fieldnames(obj);
                  tmp = strnearest(fieldx,list);
                  nC = numel(tmp);
                  choice = cell(1,nC);
                  for n=1:nC
                      choice{n} = cell2mat(list(tmp(n)));
                  end
                  fprintf(['Sorry, the field ',fieldx,' does not exist, you may be want one of those: ',strjoin(choice,' or '),'\n']);
                  flagerror = true;
              end
            
            if ~isprop(obj,fieldy)
                list = fieldnames(obj);
                  tmp = strnearest(fieldy,list);
                  nC = numel(tmp);
                  choice = cell(1,nC);
                  for n=1:nC
                      choice{n} = cell2mat(list(tmp(n)));
                  end
                  fprintf(['Sorry, the field ',fieldy,' does not exist, you may be want one of those: ',strjoin(choice,' or '),'\n']);
                  flagerror = true;
            end
                
            if flagerror
                return;
            end
            
            %% Get values and check if they're not empty
            F1 = getfield(obj,fieldx);
            F2 = getfield(obj,fieldy);
            
            if isempty(F1)
                fprintf(['Sorry, the field ',fieldx,' is empty\n']);
            end
            if isempty(F2)
                fprintf(['Sorry, the field ',fieldy,' is empty\n']);
            end
            
            if numel(F2)~=numel(F1)
                fprintf(['There''re not as many, ' fieldy,' values than ' fieldx,'\n']);
                flagerror = true;
                return;
            end
            
            % Check if it exist a field containing the precision of the
            % required parameters to plot
            if isprop(obj,[fieldx(1:3),'D',fieldx(4:end)])
                dF1 = getfield(obj,[fieldx(1:3),'D',fieldx(4:end)]);
            else
                dF1 = 0*F1;
            end
            if isprop(obj,[fieldy(1:3),'D',fieldy(4:end)])
                dF2 = getfield(obj,[fieldy(1:3),'D',fieldy(4:end)]);
            else
                dF2 = 0*F2;
            end

	   
            %% Apply threshold
            idx_good = 1:numel(F1);
            if ~isempty(thresx)
                idx_good = find(F1>thresx(1) & F1<thresx(2));
                F1 = F1(idx_good);
                dF1 = dF1(idx_good);
                F2 = F2(idx_good);
                dF2 = dF2(idx_good);
            end
            if ~isempty(thresy)
                idx_good = find(F2>thresy(1) & F2<thresy(2));
                F1 = F1(idx_good);
                F2 = F2(idx_good);
                dF1 = dF1(idx_good);
                dF2 = dF2(idx_good);
            end
            
            %% Retrieve the third part the user want to sort data
            nSort = 1;colors_p = mkc; idSort = {1:numel(idx_good)};
            if ~isempty(sortby)
                F3 = getfield(obj,sortby);
                if isempty(F3)
                    fprintf(['Sorry, the field ',sortby,' is empty, I can''t sort results regarding this field\n']);
                else
                    %Get and sort values
                    F3  = F3(idx_good);
                    if isnumeric(F3) && ~strcmp(sortby,'SYSGAIN')
                        F3 = roundn(F3,0);
                    end
                    vSort = unique(F3);
                    nSort = numel(vSort);                                       
                    idSort = cell(1,nSort);
                    for k=1:nSort
                        if isnumeric(F3)
                            idSort{k} = find(F3 == vSort(k));
                        else
                            idSort{k} = find(contains(F3,vSort{k}));
                        end
                    end
                end
            end
            
        
            %% Plot
            
            % Define the markers list
            if isempty(mk)
                mk = {'o','s', 'd', 'p', 'v','*', '^','>','<', '.'};
                nmk = numel(mk);
            else
                mk = {repmat(mk,[nSort,1])};
            end
            
            % Define the color scale
            if nSort <10
                o  =[1,0.5,0]; %orange
                g  =[0,0.75,0.25]; %orange
                p = [0.75,0,0.5];%purple
                y = [1,0.75,0.1]; %yellow mustard
                bb = [0,0.75,1];%blue
                colors_p = {'b','r','k','m',o,g,p,y,bb};
            else
                c1 = [0.5,0,1];
                c2 = [1,0.5,0];
                r  = linspace(c1(1),c2(1),nSort)';
                g =linspace(c1(2),c2(2),nSort)';
                b = linspace(c1(3),c2(3),nSort)';
                colors_p = [r,g,b];
            end
            

            % Plot with error bars in linear scale
            figure;
            if nSort > 1
                for k=1:nSort
                    errorbar(F1(idSort{k}),F2(idSort{k}),-dF2(idSort{k}),dF2(idSort{k}),-dF1(idSort{k}),dF1(idSort{k}),'LineStyle','none','Color',colors_p{k},'Marker',mk{max(mod(k,nmk+1),1)},'MarkerFaceColor',colors_p{k},'MarkerSize',mks);
                    hold on;
                end
            else
                errorbar(F1,F2,-dF2,dF2,-dF1,dF1,'LineStyle','none','Color',mkc,'Marker',mk{1},'MarkerFaceColor',mkc,'MarkerSize',mks);hold on;
            end
            
            %Manage the scale
            if strcmp(scale,'logy') || strcmp(scale,'loglog')
                set(gca,'YScale','log');
            elseif strcmp(scale,'logx') || strcmp(scale,'loglog')
                set(gca,'XScale','log');
            end 
            
            % Draw the line x=y
            if xyline
                xx = [xlim()];
                yy = [ylim()];
                mn = min(min([xx(:),yy(:)]));
                mx = max(max([xx(:),yy(:)]));
                plot([mn,mx],[mn,mx],'k--');
            end
            
            % Perform the polyfit
            if any(polyfitOrder)
                nFit = numel(polyfitOrder);
                coefsFit = cell(1,nFit);
                xx = linspace(min(F1(:)),max(F1(:)),1e3);
                for j=1:nFit
                    for k=1:nSort
                        xData = F1(idSort{k}).^(polyfitOrder(j));
                        yData = F2(idSort{k});
                        coefsFit{j} = polyfit(xData,yData,1);
                        plot(xx,polyval(coefsFit{j},xx.^(polyfitOrder(j))) ,'Color','k','LineStyle','-.','LineWidth',2);
                        coefsFit{j}
                    end
                end     
                
            end
            ylim([min(F2(:)),max(F2(:))]);
            % Put legends
            ylabel(legendy,'interpreter','latex','FontSize',fontsize);
            xlabel(legendx,'interpreter','latex','FontSize',fontsize);
            set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex');
            pbaspect([1.6,1,1]);
            if nSort >1
                if isnumeric(vSort)
                    legend(sprintfc('%.2f',vSort),'interpreter','latex','FontSize',fontsize,'Location',legendLoc);
                else
                    legend(vSort,'interpreter','latex','FontSize',fontsize,'Location',legendLoc);
                end
            end

            
        end
        
        function out = displayHistogram(obj,field,varargin)
            inputs = inputParser;
            inputs.addRequired('obj', @(x) isa(x,'statisticalAnalysis'));
            inputs.addRequired('field', @ischar );
            inputs.addParameter('field2', '',@ischar );
            inputs.addParameter('sortby',[], @ischar );
            inputs.addParameter('resultsof','psfr', @ischar );
            inputs.addParameter('fontsize',20, @isnumeric );
            inputs.addParameter('xlab',field, @ischar );
            inputs.addParameter('legend1','', @ischar );
            inputs.addParameter('legend2','', @ischar );
            inputs.addParameter('Color','b', @(x) isa(x,'char') || isa(x,'numeric') );
            inputs.addParameter('thres',[], @isnumeric );
            inputs.addParameter('legendLoc','northeast', @ischar );
            inputs.parse(obj,field,varargin{:});            
                                   
            %1\ Parsing inputs
            resultsof = inputs.Results.resultsof;
            fontsize = inputs.Results.fontsize;
            xlab = inputs.Results.xlab;
            mkc = inputs.Results.Color;
            thres = inputs.Results.thres;
            sortby = inputs.Results.sortby;
            field2 = inputs.Results.field2;
            legendLoc = inputs.Results.legendLoc;
            legend1 = inputs.Results.legend1;
            legend2 = inputs.Results.legend2;
            if isempty(legend1)
                legend1 = field;
            end
            if isempty(legend2)
                legend2 = field2;
            end

            if ~strcmpi(resultsof,'PSFR') && ~strcmpi(resultsof,'PRIME') && ~strcmpi(resultsof,'ALL')
                fprintf('Sorry, I can only plot PSF-R, PRIME or both results\n');
                return
            end
            
            %2\ Check the field
            if ~isprop(obj,field)
                list = fieldnames(obj);
                tmp = strnearest(field,list);
                nC = numel(tmp);
                choice = cell(1,nC);
                for n=1:nC
                    choice{n} = cell2mat(list(tmp(n)));
                end
                fprintf(['Sorry, the field ',field,' does not exist, you may be want one of those: ',strjoin(choice,' or '),'\n']);
                return
            end

            if ~isempty(field2) && ~isprop(obj,field2)
                list = fieldnames(obj);
                tmp = strnearest(field,list);
                nC = numel(tmp);
                choice = cell(1,nC);
                for n=1:nC
                    choice{n} = cell2mat(list(tmp(n)));
                end
                fprintf(['Sorry, the field ',field,' does not exist, you may be want one of those: ',strjoin(choice,' or '),'\n']);
                return
            end

            %3\ Get values and check if they're not empty
            F1 = getfield(obj,field);
            if isempty(F1)
                fprintf(['Sorry, the field ',field,' is empty\n']);
                return
            end
            
            if ~isempty(field2)
                F2 =  getfield(obj,field2);
                if isempty(F2)
                    fprintf(['Sorry, the field ',field,' is empty\n']);
                    return
                end
            end
            %4\ Apply threshold
            idx_good = 1:numel(F1);
            if ~isempty(thres)
                idx_good = find(F1>thres(1) & F1<thres(2));
                F1 = F1(idx_good);
            end
            if ~isempty(field2)
                F2 = F2(idx_good);
                idx_good = find(F2>thres(1) & F2<thres(2));
                F2 = F2(idx_good);
                F1 = F1(idx_good);
            end

            %4\ Retrieve the third part the user want to sort data
            nSort = 1;colors_p = mkc; idSort = {1:numel(idx_good)};
            if ~isempty(sortby)
                F3 = getfield(obj,sortby);
                if isempty(F3)
                    fprintf(['Sorry, the field ',sortby,' is empty, I can''t sort results regarding this field\n']);
                else
                    %Get and sort values
                    F3  = F3(idx_good);
                    if isnumeric(F3) && ~strcmp(sortby,'SYSGAIN')
                        F3 = roundn(F3,0);
                    end
                    vSort = unique(F3);
                    nSort = numel(vSort);                 
                    idSort = cell(1,nSort);
                    for k=1:nSort
                        if isnumeric(F3)
                            idSort{k} = find(F3 == vSort(k));
                        else
                            idSort{k} = find(contains(F3,vSort{k}));
                        end
                    end
                end
            end
            
            %manage colors
            if nSort <10
                o  =[1,0.5,0]; %orange
                g  =[0,0.75,0.25]; %orange
                p = [0.75,0,0.5];%purple
                y = [1,0.75,0.1]; %yellow mustard
                bb = [0,0.75,1];%blue
                colors_p = {'b','r','k','m',o,g,p,y,bb};
            else
                c1 = [0.5,0,1];
                c2 = [1,0.5,0];
                r  = linspace(c1(1),c2(1),nSort)';
                g =linspace(c1(2),c2(2),nSort)';
                b = linspace(c1(3),c2(3),nSort)';
                colors_p = [r,g,b];
            end

            %5\ Display
            nX = sqrt(nSort);
            nX = floor(nX);
            nY = ceil(nSort/nX);
            figure;
            if nSort > 1
                for k=1:nSort
                    subplot(nX,nY,k)
                    histogram(F1(idSort{k}),'FaceColor',colors_p{1});
                    if ~isempty(field2)
                        hold on;
                        histogram(F2(idSort{k}),'FaceColor',colors_p{2});          
                        legend({legend1,legend2},'interpreter','latex','fontsize',fontsize,'Location',legendLoc)
                    end
                    ylabel('Counts','interpreter','latex','FontSize',fontsize);
                    xlabel(xlab,'interpreter','latex','FontSize',fontsize);
                    set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex');
                    pbaspect([1,1,1]);
                    ymax = ylim();
                    ylim(ymax*1.1);
                    if isnumeric(vSort)
                        title(sprintfc('%d',vSort(k)),'interpreter','latex','FontSize',fontsize);
                    else
                        title(vSort{k},'interpreter','latex','FontSize',fontsize);
                    end
                end
            else
                histogram(F1,'FaceColor',mkc);
                if ~isempty(field2)
                    hold on;
                    histogram(F2);
                    legend({legend1,legend2},'interpreter','latex','fontsize',fontsize,'Location',legendLoc)
                end
                ylabel('Counts','interpreter','latex','FontSize',fontsize);
                xlabel(xlab,'interpreter','latex','FontSize',fontsize);
                set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex');
                pbaspect([1.6,1,1]);
                ymax = ylim();
                ylim(ymax*1.1);
            end
            
           if nargout
               if nSort > 1
                   out = zeros(2,nSort+1);
                   
                   for k=1:nSort
                       tmp = F1(idSort{k});
                       out(1,k) = median(tmp(:));
                       out(2,k)  = std(tmp(:));
                   end
                   out(1,end) = median(F1(:));
                   out(2,end) = std(F1(:));
               else
                   out=  [median(F1(:)) std(F1(:))];
               end
                   
               
           end
            
                      
        end
        
        function out = displayAOerrorBreakdown(obj,varargin)
            inputs = inputParser;
            inputs.addRequired('obj', @(x) isa(x,'statisticalAnalysis'));
            inputs.addParameter('obj_name',[],@ischar );
            inputs.addParameter('date',[],@ischar );
            inputs.addParameter('resultsof','psfr',@ischar );
            inputs.addParameter('sortby',[],@ischar );
            inputs.addParameter('fontsize',20,@isnumeric );
            inputs.addParameter('separatedFig',true,@islogical );
            inputs.parse(obj,varargin{:});
            
            resultsof = inputs.Results.resultsof;
            sortby = inputs.Results.sortby;
            fontsize = inputs.Results.fontsize;
            separatedFig = inputs.Results.separatedFig;
            
            % Single case ?
            objname = inputs.Results.obj_name;
            date = inputs.Results.date;
            if ~isempty(objname) && ~isempty(date)
                %1\ Check the object name
                idx = (strcmp(obj.obj_name,objname));
                if ~any(idx)
                    fprintf('The object does not exist\n');
                    return
                end
                %2\ Check the date
                idx = (strcmp(obj.DATE,date));
                if ~any(idx)
                    fprintf('I don''t have any date for this observing date\n');
                    return
                end
                %3\ Check if both exist simultaneously
                idx = (strcmp(obj.DATE,date) & strcmp(obj.obj_name,objname) );
                if ~any(idx)
                    fprintf('This file has not been acquired during this observing night'\n');
                    return
                end
                %4\ Grab the data
                id = find(idx);  
                if strcmp(resultsof,'psfr')
                    ncpa = obj.RECNCPA(id);
                    fit = obj.RECFIT(id);
                    lag = obj.RECLAG(id);
                    noise = obj.RECNOI(id);
                    alias = obj.RECALIAS(id);
                    aniso = obj.RECANISO(id);
                    tt = obj.RECTT(id);
                    noiseTT = obj.RECNOITT(id);
                else
                    ncpa = obj.PRINCPA(id);
                    fit = obj.PRIFIT(id);
                    lag = obj.PRILAG(id);
                    noise = obj.PRINOI(id);
                    alias = obj.PRIALIAS(id);
                    aniso = obj.PRIANISO(id);
                    tt = obj.PRITT(id);
                    noiseTT = obj.PRINOITT(id);                    
                end                                                               
                out={ncpa,fit,lag,noise,alias,aniso,tt,noiseTT};
                figure;
                statisticalAnalysis.displayPieChart(out,'fontsize',fontsize);
                title([objname,' - ',date,newline],'interpreter','latex','fontsize',fontsize+8,'color','b');
            else
            
                if ~isempty(sortby)
                    F3 = getfield(obj,sortby);
                    if isempty(F3)
                        fprintf(['Sorry, the field ',sortby,' is empty, I can''t sort results regarding this field\n']);
                        return;
                    else
                        %Get and sort values
                        if isnumeric(F3) && ~strcmp(sortby,'SYSGAIN')
                            F3 = roundn(F3,0);
                        end
                        vSort = unique(F3);
                        nSort = numel(vSort);
                        idSort = cell(1,nSort);
                        for k=1:nSort
                            if isnumeric(F3)
                                idSort{k} = find(F3 == vSort(k));
                            else
                                idSort{k} = find(contains(F3,vSort{k}));
                            end
                        end
                    end
                    nX = sqrt(nSort);
                    nX = floor(nX);
                    nY = ceil(nSort/nX);
                else
                    idSort = {1:numel(obj.RECSRMAR)};
                    nSort = 1;
                    nX = 1;nY =1;
                end
            
                figure;
                for k=1:nSort
                    if strcmp(resultsof,'psfr')                       
                        ncpa = obj.RECNCPA(idSort{k});
                        fit = obj.RECFIT(idSort{k});
                        lag = obj.RECLAG(idSort{k});
                        noise = obj.RECNOI(idSort{k});
                        alias = obj.RECALIAS(idSort{k});
                        aniso = obj.RECANISO(idSort{k});
                        tt = obj.RECTT(idSort{k});
                        noiseTT = obj.RECNOITT(idSort{k});
                    else                        
                        ncpa = obj.PRINCPA(idSort{k});
                        fit = obj.PRIFIT(idSort{k});
                        lag = obj.PRILAG(idSort{k});
                        noise = obj.PRINOI(idSort{k});
                        alias = obj.PRIALIAS(idSort{k});
                        aniso = obj.PRIANISO(idSort{k});
                        tt = obj.PRITT(idSort{k});
                        noiseTT = obj.PRINOITT(idSort{k});
                    end
                    
                    out={ncpa,fit,lag,noise,alias,aniso,tt,noiseTT};
                    
                    %\5 Display the pie chart
                    if ~separatedFig
                        subplot(nX,nY,k)
                    end
                    statisticalAnalysis.displayPieChart(out,'fontsize',fontsize);
                    if ~isempty(sortby)
                        if isnumeric(vSort)
                            title([sprintfc('%.2f',vSort(k)),newline],'interpreter','latex','FontSize',fontsize+8,'color','b');
                        else
                            title([vSort(k),newline],'interpreter','latex','fontsize',fontsize+8,'color','b');
                        end                                        
                    end
                    
                    if separatedFig && k<nSort
                        figure;
                    end                 
                end
            end
        end
        
        
        
    end
    
    methods (Static)
        
        function displayPieChart(val,varargin)
            inputs = inputParser;
            inputs.addRequired('val', @iscell);
            inputs.addParameter('sortby',[],@ischar );
            inputs.addParameter('fontsize',16,@isnumeric );
            inputs.parse(val,varargin{:});
            fontsize = inputs.Results.fontsize;
            
            % Grab Median values     
            ncpa = mean(val{1}); fit = mean(val{2});
            lag = mean(val{3}); noise = mean(val{4});
            alias = mean(val{5}); aniso = mean(val{6});
            tt = mean(val{7}); noiseTT = mean(val{8});
            
            % Grab std values
            nC = length(val{1});
            if nC > 1
                ncpastd = std(val{1}); fitstd = std(val{2});
                lagstd = std(val{3}); noisestd = std(hypot(val{4},val{8}));
                aliasstd = std(val{5}); anisostd = std(val{6});
                ttstd = std(val{7}); 
            end
            
            
            % Define values and legend vectors
            vec = [ncpa, fit,lag,tt,hypot(noise,noiseTT),alias,aniso];
            nV = numel(vec);
            if nC > 1
                leg = {...
                    ['NCPA',newline,num2str(ncpa,3),' nm'],...
                    ['Fitting',newline,num2str(fit,3),' nm',newline '$\sigma$=' num2str(fitstd,3),' nm'] ,...
                    ['Lag',newline,num2str(lag,3),' nm',newline '$\sigma$=' num2str(lagstd,3),' nm'] , ...
                    ['Jitter',newline,num2str(tt,3),' nm',newline '$\sigma$=' num2str(ttstd,3),' nm'],...
                    ['Noise',newline,num2str(hypot(noise,noiseTT),3),' nm',newline '$\sigma$=' num2str(noisestd,3),' nm'],...
                    ['Aliasing',newline,num2str(alias,3), 'nm',newline '$\sigma$=' num2str(aliasstd,3),' nm'],...
                    ['Aniso',newline,num2str(aniso,3), 'nm',newline '$\sigma$=' num2str(anisostd,3),' nm']};
            else
                leg = {...
                    ['NCPA',newline,num2str(ncpa,3),' nm'],...
                    ['Fitting',newline,num2str(fit,3),' nm'] ,...
                    ['Lag',newline,num2str(lag,3),' nm'] , ...
                    ['Jitter',newline,num2str(tt,3),' nm'],...
                    ['Noise',newline,num2str(hypot(noise,noiseTT),3),' nm'],...
                    ['Aliasing',newline,num2str(alias,3), 'nm'],...
                    ['Aniso',newline,num2str(aniso,3), 'nm']};
            end
            
            %PIE chart
            pp = pie(vec,false(1,nV),leg);
            for k=1:nnz(vec)
                pp(2*k).Interpreter = 'latex';
                pp(2*k).FontSize = fontsize;
            end
            
        end
    end
end
