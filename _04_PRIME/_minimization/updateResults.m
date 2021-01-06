function obj = updateResults(obj)
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'prime'));     
            inputs.parse(obj);
            
            %1\ Grab results
            beta_   = obj.beta_;
            J       = obj.JacobianMatrix_;
            fRes    = obj.residualMap_;
            nStars  = obj.psfr.trs.src(obj.idSrc).nObj;
            
            %2\ Compensate photometry for the normalization factor
            beta_(1:nStars) = beta_(1:nStars).*obj.normFactor;
            
            %3\ Get precision of estimates
            dx = diff(nlparci(real(beta_),fRes,'jacobian',J),1,2);
            if any(dx~=dx)
                for i=1:length(dx)
                    dx(i) = diff(nlparci(real(beta_(i)),fRes,'jacobian',J(:,i)),1,2);
                end
            end
            dx              = reshape(dx,1,[]);
            obj.x_final     = beta_;
            obj.x_prec      = dx;
            
            %4\ Unpacking atmosphere parameters     
            d                       = nStars*3;
            obj.atm_fit             = [];            
            obj.atm_fit.wavelength  = obj.psfr.trs.cam.wavelength;
           
            % r0 estimation
            idx                     = obj.idxR0 +d;
            obj.atm_fit.r0          = sum(beta_(idx))^(-3/5);
            obj.atm_fit.dr0         = 3/5*sum(beta_(idx))^(-8/5)*sqrt(sum(dx(idx).^2));
            obj.atm_fit.r0_500nm    = obj.atm_fit.r0*(0.5e-6/obj.psfr.trs.cam.wavelength)^1.2;
            % Cn2 estimation
            idx                     = obj.idxCn2 +d;
            if isempty(obj.idxCn2)
                obj.atm_fit.nLayer  = 1;
                obj.atm_fit.weights = 1;
                obj.atm_fit.heights = 0;               
                obj.atm_fit.Cn2     = obj.atm_fit.r0;
                obj.atm_fit.dCn2    = dx(obj.idxR0 +d);
            else
                obj.atm_fit.nLayer  = numel(beta_(idx));
                obj.atm_fit.weights = beta_(idx)/sum(beta_(idx));
                obj.atm_fit.heights = [obj.psfr.trs.atm.heights];
                obj.atm_fit.Cn2     = beta_(idx);
                obj.atm_fit.dCn2    = dx(idx);
            end
            
                       
            % 5\ Unpacking gains parameters
            obj.gains_fit.gAO   = obj.x_fixed.gAO ;
            obj.gains_fit.dgAO  = 0;
            obj.gains_fit.gTT   = obj.x_fixed.gTT;
            obj.gains_fit.dgTT  = 0;
            obj.gains_fit.gAl   = obj.x_fixed.gAl;
            obj.gains_fit.dgAl  = 0;
            if ~isempty(obj.idxDao)
                obj.gains_fit.gAO   = beta_(obj.idxDao+d);
                obj.gains_fit.dgAO  = dx(obj.idxDao+d);
            end
            if ~isempty(obj.idxDtt)
                obj.gains_fit.gTT   = beta_(obj.idxDtt+d);
                obj.gains_fit.dgTT  = dx(obj.idxDtt+d);
            end
            if ~isempty(obj.idxDal)
                obj.gains_fit.gAl   = beta_(obj.idxDal+d);
                obj.gains_fit.dgAl  = dx(obj.idxDal+d);
            end
                       
            % 6\ Getting the fitted map
            if ~isempty(obj.idxStatModes)
                fact = 1e-9*2*pi/obj.psfr.trs.cam.wavelength;
                obj.map_fit.coefs           = beta_(obj.idxStatModes+d);
                obj.map_fit.dcoefs          = dx(obj.idxStatModes+d);
                obj.map_fit.map_zer         = fact*obj.psfr.trs.tel.pupil.*reshape(obj.statModes*obj.map_fit.coefs',obj.psfr.trs.tel.resolution,[]);
                obj.map_fit.map_stat        = obj.map_fit.map_zer + obj.psfr.res.static.stat_map_interp*fact;
                obj.map_fit.wfsFocusCalib   = 0.198*obj.psfr.trs.wfs.zstage_defocus; %micron/mm
                obj.map_fit.wfsFocusMeas    = obj.map_fit.coefs(obj.map_fit.jindex==4)*1e-3;
                obj.map_fit.dwfsFocusMeas   = obj.map_fit.dcoefs(obj.map_fit.jindex==4)*1e-3;
            end
            
            %7\ Unpacking and concatenates Stellar parameters
            xstars_fit              = beta_(1:d);
            xstars_prec             = dx(1:d);
            obj.catalog_fit.id      = 1:nStars;
            obj.catalog_fit.x       = xstars_fit(nStars+1:nStars*2)*obj.psfr.trs.cam.pixelScale;
            obj.catalog_fit.y       = xstars_fit(2*nStars+1:nStars*3)*obj.psfr.trs.cam.pixelScale;
            obj.catalog_fit.flux    = xstars_fit(1:nStars);
            obj.catalog_fit.dx      = xstars_prec(nStars+1:nStars*2)*obj.psfr.trs.cam.pixelScale;
            obj.catalog_fit.dy      = xstars_prec(2*nStars+1:nStars*3)*obj.psfr.trs.cam.pixelScale;
            obj.catalog_fit.dflux   = xstars_prec(1:nStars);
            obj.catalog_fit.dmag    = 2.5*obj.catalog_fit.dflux./obj.catalog_fit.flux/log(10);
            if obj.fitBg
                beta_(end)          = beta_(end).*obj.normFactor;
                obj.catalog_fit.bg  = beta_(end);
                obj.catalog_fit.dbg = dx(end).*obj.normFactor;
            end
            
            
            %8\ Grab the best-fitted PSF
            wMap            = obj.weightMap;
            obj.weightMap   = 1;  % Put the weighting map to ones to have the reconstruction over all pixels          
            obj.rec_        = imageModel(beta_,obj.xdata,obj);
            obj.psf         = psfStats(obj.rec_,obj.psfr.trs.tel.pupil,obj.psfr.trs.cam.wavelength,obj.psfr.trs.cam.samp,...
                obj.psfr.trs.cam.pixelScale,obj.psfr.trs.src(obj.idSrc));
            obj.weightMap   = wMap;              
            
            %9\ Get the residual error
            obj.catalog_fit.fvu = puakoTools.getFVU(obj.psfr.trs.cam.image(:,:,obj.idSrc),obj.psf.image);
        end
        
        