function obj = fittingSetup(obj,varargin)
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'prime'));
            inputs.addParameter('aoinit',[],@isnumeric);
            inputs.addParameter('aobounds',[],@isnumeric);
            inputs.addParameter('fitR0',true,@islogical);
            inputs.addParameter('fitCn2',false,@islogical);
            inputs.addParameter('fitGains',[true,true,true],@islogical);
            inputs.addParameter('nZernMode',[],@isnumeric);
            inputs.parse(obj,varargin{:});
            
            nZernMode_history = obj.psfr.trs.mat.nZernMode;
            obj.nZernMode= inputs.Results.nZernMode;
            fitCn2       = inputs.Results.fitCn2;
            fitR0        = inputs.Results.fitR0;
            fitGains     = inputs.Results.fitGains;
            aoinit_      = inputs.Results.aoinit;
            aobounds_    = inputs.Results.aobounds;
            
            
            %1\ Recalculate the Residual phase structure function if the
            %user opts for a modal gain retrieval            
            if isempty(obj.nZernMode)
                obj.psfr.sf.Dho_z = obj.psfr.sf.Dho;
                obj.psfr.cov.Cho_z = obj.psfr.cov.Cho;
                obj.nGainsHO = 1;
            elseif  (~isempty(obj.nZernMode)) && (numel(obj.nZernMode) ~= numel(nZernMode_history) || (any(obj.nZernMode ~= nZernMode_history)))
                [obj.psfr.sf.Dho_z,obj.psfr.cov.Cho_z] = computeResidualPhaseStructureFunction(obj.psfr,'nZernMode',obj.nZernMode);
                obj.nGainsHO = numel(obj.nZernMode)+1; % if 1 the gain is common to each mode
            end
            
            %2\ Define initial guesses
            if (~isempty(aoinit_)) && (~isempty(aobounds_))
                ao_init = aoinit_;
                ao_lb   = aobounds_(1,:);
                ao_ub   = aobounds_(2,:);
            else                                
                obj.idxCn2 = [];Cn2_init  = []; Cn2_lb    = [];
                Cn2_ub    = [];obj.idxR0  = [];                
                gain_init = [];gain_lb   = [];gain_ub   = [];                                
                nInit    = 1;       
                
                %2.1 Cn2 and r0 init                              
                if fitCn2
                    obj.idxCn2= nInit:nInit+obj.atm.nLayer-1;
                    obj.idxR0 = obj.idxCn2;
                    nInit     = nInit+obj.atm.nLayer;
                    Cn2_init  = [obj.atm.layer.fractionnalR0]*obj.atm.r0^(-5/3)*(obj.atm.wavelength/obj.wvl_ref)^2;
                    Cn2_lb    = 5^(-5/3)*ones(1,obj.atm.nLayer);
                    Cn2_ub    = 0.05^(-5/3)*ones(1,obj.atm.nLayer);
                else
                    if fitR0
                        obj.idxR0 = nInit;
                        nInit     = nInit+1;
                        Cn2_init  = obj.psfr.trs.atm.r0^(-5/3)*(obj.psfr.trs.atm.wavelength/obj.psfr.trs.cam.wavelength)^2;
                        Cn2_lb    = 5^(-5/3);
                        Cn2_ub    = 0.01^(-5/3);
                    end
                end
                
                %2.2. HODM gain              
                if fitGains(1)
                    obj.idxDho = nInit:nInit+numel(obj.nZernMode);
                    nInit          = nInit+length(obj.idxDho);
                    gain_init    = ones(1,numel(obj.idxDho));
                    gain_lb      = [gain_lb,zeros(1,numel(obj.idxDho))];
                    gain_ub     = [gain_ub,1e2*ones(1,numel(obj.idxDho))];               
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
                    gain_init    = [gain_init,1];
                    gain_lb      = [gain_lb,0];
                    gain_ub      = [gain_ub,1e2];                
                end
                
                %2.5 Concatenating vectors
                ao_init   = [Cn2_init,gain_init];
                ao_lb     = [Cn2_lb,gain_lb];
                ao_ub     = [Cn2_ub,gain_ub];
            end
            
            %3\  Define initial guess on stellar parameters
            uMax = 2; % PSF position may change by +/- 2 pixels
            photo_init = ones(1,obj.nStars)/obj.nStars;            
            pos_init   = [0,0];%1e3*[obj.psfr.trs.src.x,obj.psfr.trs.src.y]/obj.psfr.trs.cam.pixelScale;
            bg_init    = 0;
            photo_b   = [zeros(1,obj.nStars);2*ones(1,obj.nStars)];
            pos_b      = uMax*ones(1,2*obj.nStars);
            bg_b       = 5*std(reshape(obj.im_sky,obj.fov_sky^2,[]),1);
            
            %4\ Concatenating PSF and stellar parameters
            obj.x_init = [photo_init,pos_init,ao_init,bg_init];
            obj.lbounds= [photo_b(1,:),pos_init-pos_b,ao_lb,-bg_b];
            obj.ubounds= [photo_b(2,:),pos_init+pos_b,ao_ub,bg_b];
            
            %5\ Define the additional useful data
            obj.xdata = {obj.idxCn2,obj.idxR0,obj.idxDho,obj.idxDtt,obj.idxDal};
        end