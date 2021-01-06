function obj = fittingSetup(obj,varargin)
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'prime'));
            inputs.addParameter('aoinit',[],@isnumeric);
            inputs.addParameter('aobounds',[],@isnumeric);
            inputs.addParameter('fitR0',true,@islogical);
            inputs.addParameter('fitCn2',false,@islogical);
            inputs.addParameter('fitGains',[true,true,true],@islogical);
            inputs.addParameter('fitBg',false,@islogical);
            inputs.addParameter('jZernGain',[],@isnumeric);
            inputs.addParameter('fitStatModes',[],@isnumeric);
            inputs.addParameter('statModesFunction',[],@ischar);
            inputs.addParameter('umax',5,@isnumeric);
            inputs.parse(obj,varargin{:});
            
            jZernGain_history = obj.psfr.trs.mat.jZernGain;
            obj.jZernGain= inputs.Results.jZernGain;
            fitCn2       = inputs.Results.fitCn2;
            fitR0        = inputs.Results.fitR0;
            fitGains     = inputs.Results.fitGains;
            if numel(fitGains) == 1
                fitGains = [fitGains fitGains fitGains];
            end
            fitStatModes = inputs.Results.fitStatModes;
            aoinit_      = inputs.Results.aoinit;
            aobounds_    = inputs.Results.aobounds;
            fitBg        = inputs.Results.fitBg;
            statModesFunction = inputs.Results.statModesFunction;
            umax         = inputs.Results.umax;
            
            %1\ Recalculate the Residual phase structure function if the
            %user opts for a modal gain retrieval            
            if isempty(obj.jZernGain)
                obj.psfr.sf.Dao_z   = obj.psfr.sf.Dao;
                obj.psfr.cov.Cao_z  = obj.psfr.cov.Cao;
                obj.nGainsHO = 1;
            elseif  (~isempty(obj.jZernGain)) && (numel(obj.jZernGain) ~= numel(jZernGain_history) || (any(obj.jZernGain ~= jZernGain_history)))
                obj.psfr        = computeResidualPhaseStructureFunction(obj.psfr,'jZernGain',obj.jZernGain);
                obj.nGainsHO    = numel(obj.jZernGain)+1; % if 1 the gain is common to each mode
            end
            
            %2\ Define initial guesses
            if (~isempty(aoinit_)) && (~isempty(aobounds_))
                ao_init = aoinit_;
                ao_lb   = aobounds_(1,:);
                ao_ub   = aobounds_(2,:);
            else                                
                obj.idxCn2  = [];Cn2_init   = []; Cn2_lb    = [];
                Cn2_ub      = [];obj.idxR0  = [];                
                gain_init   = [];gain_lb    = []; gain_ub   = [];        
                statCoefs_init = [];statCoefs_lb = []; statCoefs_ub = [];
                nInit       = 1;       
                
                %2.1 Cn2 and r0 init           
                r053 = (obj.psfr.trs.res.seeing.r0*(obj.psfr.trs.cam.wavelength/0.5e-6)^1.2)^(-5/3);
                if fitCn2
                    nL = obj.psfr.trs.atm.nLayer;
                    obj.idxCn2= nInit:nInit+nL-1;
                    obj.idxR0 = obj.idxCn2;
                    nInit     = nInit+nL;
                    Cn2_init  = obj.psfr.trs.atm.weights*r053;
                    Cn2_lb    = zeros(1,nL);
                    Cn2_ub    = 0.05^(-5/3)*ones(1,nL);
                else
                    if fitR0
                        obj.idxR0 = nInit;
                        nInit     = nInit+1;
                        Cn2_init  = r053;
                        Cn2_lb    = 0;
                        Cn2_ub    = 0.01^(-5/3);
                    end
                end
                
                %2.2. HODM gain              
                if fitGains(1)
                    obj.idxDao  = nInit:nInit+numel(obj.jZernGain);
                    nInit       = nInit+length(obj.idxDao);
                    gain_init   = ones(1,numel(obj.idxDao));
                    gain_lb     = [gain_lb,zeros(1,numel(obj.idxDao))];
                    gain_ub     = [gain_ub,1e2*ones(1,numel(obj.idxDao))];               
                end                
                
                %2.3. TTDM gain
                if fitGains(2)
                    obj.idxDtt   = nInit;
                    nInit        = nInit+1;
                    gain_init    = [gain_init,1];
                    gain_lb      = [gain_lb,0];
                    gain_ub      = [gain_ub,1e2];              
                end
                
                %2.4. Aliasing model gain
                if fitGains(3)
                    obj.idxDal = nInit;
                    nInit        = nInit+1;
                    gain_init    = [gain_init,1];
                    gain_lb      = [gain_lb,0];
                    gain_ub      = [gain_ub,1e2];                
                end
                
                %2.5. Static modes
                if (~isempty(fitStatModes)   ||  (numel(obj.list_fixed) > 1 && contains(upper(obj.list_fixed{1:2:end}),'AZ'))) && ~isempty(statModesFunction)
                    if isempty(fitStatModes)
                        errordlg('Please, defined the field fitStatModes as a vector of valid indexes');                
                        obj.flagError = 1;
                        return;
                    end
                    nModes = numel(fitStatModes);
                    
                    if strcmpi(statModesFunction,'ZERNIKE')
                        obj.statModes  = zernike_puako(fitStatModes,obj.psfr.trs.tel.resolution);
                        obj.statModes  = obj.statModes.modes;
                    elseif strcmpi(statModesFunction,'PETAL')
                        tmp = keckPetalModes('resolution',obj.psfr.trs.tel.resolution);
                        obj.statModes = tmp.modes(:,fitStatModes);
                    elseif strcmpi(statModesFunction,'PISTON')
                        tmp = fitsread('keckPistonModes_200x200.fits');
                        tmp = puakoTools.interpolate(tmp,obj.psfr.trs.tel.resolution,'nearest');
                        tmp = reshape(tmp,obj.psfr.trs.tel.resolution^2,36);
                        %tmp = keckPetalModes('resolution',obj.psfr.trs.tel.resolution,'petal',false);
                        obj.statModes = tmp(:,fitStatModes);
                        % manage rotation
                        for k=1:numel(fitStatModes)
                            tmp = puakoTools.rotateIm(reshape(obj.statModes(:,k),obj.psfr.trs.tel.resolution,[]),obj.psfr.trs.wfs.pupilAngle);
                            obj.statModes(:,k) = tmp(:);
                        end
                    end
                    
                    obj.map_fit.jindex  = fitStatModes;
                    statCoefs_init      = zeros(1,nModes);
                    sMax                = obj.psfr.trs.cam.wavelength*1e9;
                    statCoefs_lb        =-sMax*ones(1,nModes); %-1micron min
                    statCoefs_ub        = sMax*ones(1,nModes); %1micron max
                    
                    if (numel(obj.list_fixed) > 1 && contains(upper(obj.list_fixed{1:2:end}),'AZ'))       
                        obj.idxStatModes  = [];
                    else
                        obj.idxStatModes  = nInit:nInit+nModes-1;                                            
                    end
                end
                %2.6 Concatenating vectors
                ao_init = [Cn2_init,gain_init,statCoefs_init];
                ao_lb   = [Cn2_lb,gain_lb,statCoefs_lb];
                ao_ub   = [Cn2_ub,gain_ub,statCoefs_ub];
            end
            
            %3\  Define initial guess on stellar parameters
            nStars      = obj.psfr.trs.src.nObj;
            % flux
            photo_init  = [obj.psfr.trs.src(obj.idSrc).F]/sum(obj.psfr.trs.src(obj.idSrc).F);           
            photo_b     = [zeros(1,nStars);1e2*ones(1,nStars)];
            % positions

            %             idM         = find(obj.psfr.trs.src(obj.idSrc).F == max(obj.psfr.trs.src(obj.idSrc).F));
%             pos_init    = 1e3*[obj.psfr.trs.src(obj.idSrc).y - obj.psfr.trs.src(obj.idSrc).y(idM)...
%                                 ,obj.psfr.trs.src(obj.idSrc).x - obj.psfr.trs.src(obj.idSrc).x(idM)]...
%                                 /obj.psfr.trs.cam.pixelScale;           %in pixel            
%             [yref,xref] = find(obj.ydata == max(obj.ydata(:)));
%             pos_init    = pos_init - (obj.psfr.trs.cam.resolution/2+1 - [xref*ones(1,nStars),yref*ones(1,nStars)]);   
%             
            idMax                   = find(obj.psfr.trs.src(obj.idSrc).F == max(obj.psfr.trs.src(obj.idSrc).F));
            idMin                   = find(obj.psfr.trs.src(obj.idSrc).F ~= max(obj.psfr.trs.src(obj.idSrc).F));
            pos_init                = -[abs(obj.psfr.trs.src(obj.idSrc).y(idMin) - obj.psfr.trs.src(obj.idSrc).y(idMax)),0,...
                abs(obj.psfr.trs.src(obj.idSrc).x(idMin) - obj.psfr.trs.src(obj.idSrc).x(idMax)),0];
            pos_init                = 1e3*pos_init/obj.psfr.trs.cam.pixelScale;
            [yref,xref]             = find(obj.ydata == max(obj.ydata(:)));
            pos_init                = (pos_init - (obj.psfr.trs.cam.resolution/2+1 - [yref*ones(1,nStars),xref*ones(1,nStars)]));
            
            pos_b       = umax*ones(1,2*nStars);
            % concatenate
            obj.x_init  = [photo_init,pos_init,ao_init];
            obj.lbounds = [photo_b(1,:),pos_init-pos_b,ao_lb];
            obj.ubounds = [photo_b(2,:),pos_init+pos_b,ao_ub];
            
            if fitBg
                bg_init    = 0;
                bg_b       = 1e10;
                obj.x_init = [obj.x_init,bg_init];
                obj.lbounds= [obj.lbounds,-bg_b];
                obj.ubounds= [obj.ubounds,bg_b];
            end
                                  
            %5\ Define the additional useful data
            obj.xdata = {obj.idxCn2,obj.idxR0,obj.idxDao,obj.idxDtt,obj.idxDal,obj.idxStatModes};
            
            %6\ Define fixed parameters
           obj.x_fixed = getFixedParameters(obj.list_fixed,'Cn2',obj.psfr.trs.atm.Cn2*(0.5e-6/obj.psfr.trs.cam.wavelength)^2,'az',statCoefs_init);
            
            % recalculate the Static OTF if the user wants to work with a
            % different static map                       
            if obj.x_fixed.az~=0 ||  ~isequal(logical(obj.x_fixed.map),obj.psfr.flags.staticMap)
                obj.psfr = computeStaticOpticalTransferFunction(obj.psfr,'addStatModes',{obj.x_fixed.az,obj.statModes},'flagStaticMap',obj.x_fixed.map);
            end
            
           
        end
