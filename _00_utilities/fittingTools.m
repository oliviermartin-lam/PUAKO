classdef fittingTools < handle
    
    methods (Static)
        
        % PSF-fitting
        function [out,imFit,iminit] = findStellarParameters(im,psf,x0,varargin)
            inputs = inputParser;
            inputs.addRequired('im',@isnumeric);
            inputs.addRequired('psf',@isnumeric);
            inputs.addRequired('x0',@isnumeric);
            inputs.addParameter('ron',0,@isnumeric);
            inputs.addParameter('u_',[],@isnumeric);
            inputs.addParameter('v_',[],@isnumeric);
            inputs.addParameter('umax',5,@isnumeric);
            inputs.addParameter('MaxIter',3e2,@isnumeric);
            inputs.addParameter('TolFun',1e-18,@isnumeric);
            inputs.addParameter('MaxFunEvals',1e3,@isnumeric);
            inputs.addParameter('TolX',1e-18,@isnumeric);
            inputs.parse(im,psf,x0,varargin{:});
            ron     = inputs.Results.ron;
            u_      = inputs.Results.u_;
            v_      = inputs.Results.v_;
            umax    = inputs.Results.umax;
            MaxIter = inputs.Results.MaxIter;
            TolFun  = inputs.Results.TolFun;
            MaxFunEvals = inputs.Results.MaxFunEvals;
            TolX = inputs.Results.TolX;


            %1. DATA
            ydata       = im;
            normFactor  = sum(ydata(ydata>0));
            weightMap   = ydata~=0;
            if ron
                weightMap = weightMap./sqrt(max(ydata,0) + ron^2);
            end
            ydata = ydata.*weightMap/normFactor;
            nIm   = size(im,1);
                   
            %2. MODEL
            
            %2.1 Model definition
            model = @(x,xdata) weightMap.*puakoTools.stellarFieldModel(x,xdata,nIm,'u_',u_,'v_',v_);
            nParam = length(x0);
            out     = [];
            imFit   = [];
            beta    = [];
            %2.2 Fitting options and initial guess
            
            if mod(nParam,3)
                % maybe a background estimation
                if mod(nParam-1,3)
                    fprintf('You must provide a valid initial guess vector, which should size 3*nS or 3*nS +1 with nS the number of stars\n');
                    return
                    out   = 0;
                    imFit = 0;
                    beta  = 0;
                else
                    nS    = (length(x0)-1)/3;
                    lX    = x0(1:2*nS) - umax;
                    uX    = x0(1:2*nS) + umax;
                    lF    = zeros(1,nS);
                    uF    = 1e20*ones(1,nS);
                    lb    = [lX,lF,-1e20];
                    ub    = [uX,uF,1e20];
                end
            else
                % no background estimation
                nS    = (length(x0))/3;
                lX    = x0(1:2*nS) - umax;
                uX    = x0(1:2*nS) + umax;
                lF    = 0;
                uF    = 1e20;
                lb    = [lX,lF];
                ub    = [uX,uF];
            end
            
            opt = optimoptions(@lsqcurvefit,'MaxIter',MaxIter,...
                'TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals,...
                'InitDamping',1,'Display','iter');
            
            %3. Fitting procedure
            [beta,~,R,~,~,~,J] = lsqcurvefit(model,x0,psf,ydata,lb,ub,opt);
            
            %4. Unpacking results
            xS = beta(1:nS);
            yS = beta(nS+1:2*nS);
            fS = beta(2*nS+1:3*nS)*normFactor;
            
            %5. Measurements uncertainties
            dbeta = diff(nlparci(beta,R,'jacobian',J),1,2);
            dX    = dbeta(1:nS);
            dY    = dbeta(nS+1:2*nS);
            dF    = dbeta(2*nS+1:3*nS)*normFactor;
            
            %6. Concatenating results
        
            out = zeros(nParam,nS);
            for iS=1:nS
                out(1,iS) = xS(iS);
                out(2,iS) = yS(iS);
                out(3,iS) = fS(iS);
                out(4,iS) = dX(iS);
                out(5,iS) = dY(iS);
                out(6,iS) = dF(iS);
                if mod(nParam,3)
                    out(7,iS) = beta(end)*normFactor;
                end
            end
            
            % Fitted model
            if ~mod(nParam,3)
                res = [xS,yS,fS];
            else
                res = [xS,yS,fS,beta(end)*normFactor];
            end
            
            imFit =  puakoTools.stellarFieldModel(res,psf,nIm);
            
            if nargout > 2
                iminit = puakoTools.stellarFieldModel(x0,psf,nIm);
            end
            
        end
        
        function out = stellarFieldModel(pStars,dataModel,nIm,varargin)
            inputs = inputParser;
            inputs.addRequired('pStars',@isnumeric);
            inputs.addRequired('dataModel',@isnumeric);
            inputs.addRequired('nIm',@isnumeric);
            inputs.addParameter('u_',[],@isnumeric);
            inputs.addParameter('v_',[],@isnumeric);
            inputs.parse(pStars,dataModel,nIm,varargin{:});
            
            u_ = inputs.Results.u_;
            v_ = inputs.Results.v_;
            
            if isempty(u_)
                flagOtf = false;
            end
            
            nParam = length(pStars);
            if ~mod(nParam,3)                
                nS = nParam/3;
                xS = pStars(1:nS);
                yS = pStars(1+nS:2*nS);
                fS = pStars(1+2*nS:3*nS);
                bg = 0;
            else
                nS = (nParam-1)/3;
                xS = pStars(1:nS);
                yS = pStars(1+nS:2*nS);
                fS = pStars(1+2*nS:3*nS);
                bg = pStars(end);
            end
            
            out = 0*dataModel;
            for iS = 1:nS
                % Translate the PSF
                if flagOtf
                    % Calculate the fourier phasor
                    phasor  = exp(-2i*pi*( u_*xS(iS) + v_*yS(iS) ));
                    % Get the PSF
                    psf_i   = tools.otf2psf(dataModel.*phasor);
                else
                    psf_i = puakoTools.translateImage(dataModel,[xS(iS),yS(iS)]);
                    %psf_i = imtranslate(psfModel,[xS(iS),yS(iS)]);
                end
                % Update the image and Flux scaling:
                out = out + psf_i*fS(iS);
            end
            out = puakoTools.crop(out,nIm) + bg;
        end
        
        function out = multipleImage(x,im)
            %Grab locations
            n       = length(x)/3;
            [nx,ny] = size(im);
            xloc    = x(1:n);
            yloc    = x(n+1:2*n);
            flux    = x(2*n+1:end);
            
            if any(abs(xloc) > nx) || any(abs(yloc) > ny)
                id = abs(xloc) > nx | abs(yloc) > ny;
                nn = sum(id);
                if nn == 1
                    fprintf('Warning: 1 PSF are out the field\n');
                else
                    fprintf('Warning: %d PSFs are out the field\n',nn);
                end
                flux(id) = 0;
            end
            
            % For xloc=yloc=0 and flux=1, the procedure does modify the PSF
            % while it mustn't. To be check.
            %nx      = 2*nx;
            %ny      = 2*ny;
            out     = zeros(nx,ny);
            otf     = puakoTools.psf2otf(im);
            otf     = otf/max(otf(:));
            [u,v]   = freqspace([nx ny],'meshgrid');
            for i=1:n
                yi    = yloc(i);
                xi    = xloc(i);
                fftPhasor = exp(-1i*pi*(u*xi + v*yi));
                map_i = puakoTools.otf2psf(otf.*fftPhasor);
                % Adding images
                out = out + flux(i)*map_i/sum(map_i(:));
            end
        end
        
        % Gaussian/Moffat-fitting
        function out = gaussian(x,xdata,flagSymetric)
            if nargin < 3
                flagSymetric = false;
            end
            
            if flagSymetric
                I0 = x(1);          %Amplitude
                ax = x(2);          %x spreading
                ay = x(2);          %y-spreading
                x0 = x(3);          %x-shift
                y0 = x(4);          %y-shift
            else
                % ------- Grabbing parameters ---------%
                I0 = x(1);          %Amplitude
                ax = x(2);          %x spreading
                ay = x(3);          %y-spreading
                th = x(4)*pi/180.;  %rotation
                x0 = x(5);          %x-shift
                y0 = x(6);          %y-shift
            end
            
            %FWHM : 2*sqrt(2*log(2))*sigma
            % ------- Including shifts ---------
            X     = xdata{1};
            Y     = xdata{2};
            %Shifts
            X     = X - x0;
            Y     = Y - y0;
            %Rotation
            Xr    = X.*cos(th) + Y.*sin(th);
            Yr    = Y.*cos(th) - X.*sin(th);
            
            % Gaussian expression
            out = I0.*exp(-0.5*((Xr./ax).^2 + (Yr./ay).^2) );
        end
        
        function out = moffat(x,xdata,flagSymetric)
            
            if nargin < 3
                flagSymetric = false;
            end
            
            % ------- Grabbing parameters ---------%
            if flagSymetric
                I0 = x(1);          %Amplitude
                ax = x(2);          %x spreading
                ay = x(2);          %y-spreading
                be = x(3);          %center slope
                x0 = x(4);          %x-shift
                y0 = x(5);          %y-shift
                th = 0;
            else
                I0 = x(1);          %Amplitude
                ax = x(2);          %x spreading
                ay = x(3);          %y-spreading
                be = x(4);          %center slope
                th = x(5);          %rotation
                x0 = x(6);          %x-shift
                y0 = x(7);          %y-shift
            end
            
            % ------- Including shifts ---------
            X     = xdata{1};
            Y     = xdata{2};
            %Shifts
            X     = X - x0;
            Y     = Y - y0;
            %Rotation
            Xr    = X.*cos(th) + Y.*sin(th);
            Yr    = Y.*cos(th) - X.*sin(th);
            
            % Moffat expression
            out = I0.*(1. + (Xr./ax).^2 + (Yr./ay).^2).^(-be);
            %Normalization
            %if exist('norm','var')
            %    out = out*(be-1)/pi/ax/ay;
            %end
        end
        
        function out = multiAnalyticIsoplanatic(x,xdata,type,flagSymetric)
            if nargin < 4
                flagSymetric = false;
            end
            
            if strcmp(type,'gaussian')
                f   = @(x,xdata) puakoTools.gaussian(x,xdata,flagSymetric);
                if flagSymetric
                    nfix = 1;
                else
                    nfix = 3;
                end
                nS  = (length(x)-nfix)/3;
                xfix= x(1+nS:nfix+nS);
                fS  =  x(1:nS);
                xS  = x(nfix+1 + nS:nfix+1+2*nS);
                yS  = x(nfix+1 + 2*nS:end);
            else
                f = @(x,xdata) puakoTools.moffat(x,xdata,flagSymetric);
                if flagSymetric
                    nfix = 2;
                else
                    nfix = 4;
                end
                nS  = (length(x)-nfix)/3;
                xfix= x(1+nS:nfix+nS);
                fS  = x(1:nS);
                xS  = x(nfix+1+nS:nfix+2*nS);
                yS  = x(nfix+1+2*nS:end);
            end
            
            
            out = zeros(size(xdata,1));
            for k=1:nS
                out = out + f([fS(k) xfix,xS(k) yS(k)],xdata);
            end
        end
        
        function out = multiAnalyticAniso(x,xdata,type)
            
            if strcmp(type,'gaussian')
                f = @(x,xdata) puakoTools.gaussian(x,xdata);              
                nG  = length(x)/6;                
            else
                f = @(x,xdata) puakoTools.moffat(x,xdata);
                nG  = length(x)/7;                
            end
            
            
            out = zeros(size(xdata,1));
            for i=1:nG
                out = out + f(x(1+(i-1)*6:i*6),xdata);
            end
        end
        
        % Data fitting
        function p = polyfitweighted(x,y,n,w)
            % polyfitweighted.m
            % -----------------
            %
            % Find a least-squares fit of 1D data y(x) with an nth order
            % polynomial, weighted by w(x).
            %
            % By S.S. Rogers (2006), based on polyfit.m by The MathWorks, Inc. - see doc
            % polyfit for more details.
            %
            % Usage
            % -----
            %
            % P = polyfitweighted(X,Y,N,W) finds the coefficients of a polynomial
            % P(X) of degree N that fits the data Y best in a least-squares sense. P
            % is a row vector of length N+1 containing the polynomial coefficients in
            % descending powers, P(1)*X^N + P(2)*X^(N-1) +...+ P(N)*X + P(N+1). W is
            % a vector of weights.
            %
            % Vectors X,Y,W must be the same length.
            %
            % Class support for inputs X,Y,W:
            %    float: double, single
            %
            
            % The regression problem is formulated in matrix format as:
            %
            %    yw = V*p    or
            %
            %          3    2
            %    yw = [x w  x w  xw  w] [p3
            %                            p2
            %                            p1
            %                            p0]
            %
            % where the vector p contains the coefficients to be found.  For a
            % 7th order polynomial, matrix V would be:
            %
            % V = [w.*x.^7 w.*x.^6 w.*x.^5 w.*x.^4 w.*x.^3 w.*x.^2 w.*x w];
            
            if ~isequal(size(x),size(y),size(w))
                error('X and Y vectors must be the same size.')
            end
            
            x = x(:);
            y = y(:);
            w = w(:);
            
            
            % Construct weighted Vandermonde matrix.
            %V(:,n+1) = ones(length(x),1,class(x));
            V(:,n+1) = w;
            for j = n:-1:1
                V(:,j) = x.*V(:,j+1);
            end
            
            % Solve least squares problem.
            [Q,R] = qr(V,0);
            ws = warning('off','all');
            p = R\(Q'*(w.*y));    % Same as p = V\(w.*y);
            
            warning(ws);
            if size(R,2) > size(R,1)
                warning('polyfitweighted:PolyNotUnique', ...
                    'Polynomial is not unique; degree >= number of data points.')
            elseif condest(R) > 1.0e10
                if nargout > 2
                    warning('polyfitweighted:RepeatedPoints', ...
                        'Polynomial is badly conditioned. Remove repeated data points.')
                else
                    warning('polyfitweighted:RepeatedPointsOrRescale', ...
                        ['Polynomial is badly conditioned. Remove repeated data points\n' ...
                        '         or try centering and scaling as described in HELP POLYFIT.'])
                end
            end
            p = p.';          % Polynomial coefficients are row vectors by convention.
        end
        
    end
end