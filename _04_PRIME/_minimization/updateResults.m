function obj = updateResults(obj,beta,fRes,J)
            inputs = inputParser;
            inputs.addRequired('obj',@(x) isa(x,'prime'));
            inputs.addRequired('beta',@isnumeric);
            inputs.addRequired('fRes',@isnumeric);
            inputs.addRequired('J',@isnumeric);
            inputs.parse(obj,beta,fRes,J);
            
            
            %1\ Compensate photometry and background for the normalization factor
            beta(1:obj.nStars) = beta(1:obj.nStars).*obj.normFactor;
            beta(end) = beta(end).*obj.normFactor;
            
            %2\ Get precision of estimates
            dx = diff(nlparci(real(beta),fRes,'jacobian',J),1,2);
            for i=1:length(dx)
                dx(i) = diff(nlparci(real(beta(i)),fRes,'jacobian',J(:,i)),1,2);
            end
            dx = reshape(dx,1,[]);
            obj.x_final = beta;
            obj.x_prec  = dx;
            
            %3\ Unpacking AO parameters
            wvl_ratio =  obj.psfr.trs.cam.wavelength/500e-9;
            d = obj.nStars*3;
            % Seeing estimation
            idx         = obj.idxR0 +d;
            obj.r0_fit  = sum(beta(idx))^(-3/5);
            obj.r0_prec = 3/5*sum(beta(idx))^(-8/5)*sqrt(sum(dx(idx).^2));
            % Cn2 estimation
            idx         = obj.idxCn2 +d;
            obj.Cn2_fit = beta(idx);
            obj.Cn2_prec= dx(idx);
            if isempty(obj.Cn2_fit)
                obj.Cn2_fit = obj.r0_fit^(-5/3);
                obj.Cn2_prec= 5/3*obj.r0_prec*obj.r0_fit^(-8/3);
            end
            % Retrieved atmosphere
            if numel(obj.Cn2_fit) == 1
                hl = 0;
            else
                hl          = [obj.atm.layer.altitude];
            end
            
            nL  = length(obj.Cn2_fit*wvl_ratio^2);
            r0  = sum(obj.Cn2_fit*wvl_ratio^2)^(-3/5);
            fl  = obj.Cn2_fit/sum(obj.Cn2_fit);
           
            %bestFit.atm_fit = atmosphere(photometry.V0,r0,bestFit.atm.L0,'fractionnalR0',fl,...
              %  'altitude',hl,'layeredL0',bestFit.atm.L0,'windSpeed',wS,'windDirection',wD);
            
            % AO parameters
            obj.xao_fit = beta([obj.idxDho,obj.idxDtt,obj.idxDal]+3*obj.nStars);
            obj.xao_prec= dx([obj.idxDho,obj.idxDtt,obj.idxDal]+3*obj.nStars);
            
            %4\ Unpacking Stellar parameters
            obj.bg_fit           = beta(end);
            obj.xstars_fit       = beta(1:d);
            obj.xstars_prec      = dx(1:d);
            obj.catalog_fit.id   = 1:obj.nStars;
            obj.catalog_fit.RA   = obj.xstars_fit(obj.nStars+1:obj.nStars*2)*obj.psfr.trs.cam.pixelScale;
            obj.catalog_fit.DEC  = obj.xstars_fit(2*obj.nStars+1:obj.nStars*3)*obj.psfr.trs.cam.pixelScale;
            obj.catalog_fit.FLUX = obj.xstars_fit(1:obj.nStars);
            obj.catalog_fit.dRA  = obj.xstars_prec(obj.nStars+1:obj.nStars*2)*obj.psfr.trs.cam.pixelScale;
            obj.catalog_fit.dDEC = obj.xstars_prec(2*obj.nStars+1:obj.nStars*3)*obj.psfr.trs.cam.pixelScale;
            obj.catalog_fit.dFLUX= obj.xstars_prec(1:obj.nStars);
            obj.catalog_fit.dMAG = 2.5*obj.catalog_fit.dFLUX./obj.catalog_fit.FLUX/log(10);
            
            %5\ Grab the best-fitted PSF + 3-sigma PSF
            wMap        = obj.weightMap;
            obj.weightMap  = 1;
            obj.im.rec = imageModel(beta,obj.xdata,obj);
            if all(obj.x_prec <= max(abs(obj.lbounds),abs(obj.ubounds)))
                im_3sig_m   = obj.modelFUN(obj.x_final-obj.x_prec,obj.xdata);
                psf_3sig_m  = obj.psf_fit;
                im_3sig_p   = obj.modelFUN(obj.x_final+obj.x_prec,obj.xdata);
                psf_3sig_p  = obj.psf_fit;
                obj.im.rec_3sig = [im_3sig_m,im_3sig_p];
                obj.psf.rec_3sig= [psf_3sig_m,psf_3sig_p];
            end
            obj.weightMap= wMap;
            
        end
        
        