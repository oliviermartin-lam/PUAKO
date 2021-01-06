classdef psfReconstruction < handle
    
    properties
        % --------------------- Observation
        trs;    %Telemetry class
        % --------------------- Results
        res = [];
        otf = [];
        cov = [];
        sf = [];
        psf = [];                                                     % Reconstructed 2D PSF - psfStats class
        flags = [];
        rec_;
        statCoefs;
        statModes;
    end
    
    methods
        
        function obj = psfReconstruction(trs,varargin)
            inputs = inputParser;
            inputs.addRequired('trs',@(x) isa(x,'telemetry'));
            inputs.addParameter('fov',trs.cam.resolution*2,@isnumeric);
            inputs.addParameter('flagToeplitz',true,@islogical);
            inputs.addParameter('flagNoiseMethod','autocorrelation',@ischar);
            inputs.addParameter('flagAnisoMethod','flicker',@ischar);
            inputs.addParameter('flagAoPattern','circle',@ischar);
            inputs.addParameter('flagDphiMethod','Vii',@ischar);
            inputs.addParameter('statModes',{[]},@iscell);
            inputs.addParameter('fitbg',false,@islogical);
            inputs.parse(trs,varargin{:});
            
            %1\ Parsing inputs
            obj.trs                 = inputs.Results.trs;
            obj.otf.fov_fit         = inputs.Results.fov;
            obj.flags.toeplitz      = inputs.Results.flagToeplitz;
            obj.flags.noiseMethod   = inputs.Results.flagNoiseMethod;
            obj.flags.anisoMethod   = inputs.Results.flagAnisoMethod;
            obj.flags.aoPattern     = inputs.Results.flagAoPattern;
            obj.flags.dphiMethod    = inputs.Results.flagDphiMethod;
            obj.statModes           = inputs.Results.statModes;
            obj.flags.fitbg         = inputs.Results.fitbg;
            
            tpsfr = tic();
            %1\ Update the reconstruction prerequisites
            obj = updateReconstructionPrerequisites(obj);
            
            %2\ Reconstruct the PSF
            obj = forwardPSFR(obj,'flagToeplitz',obj.flags.toeplitz,'flagNoiseMethod',obj.flags.noiseMethod,...
                'flagAnisoMethod',obj.flags.anisoMethod,'flagAoPattern',obj.flags.aoPattern,'statModes',obj.statModes,'flagDphiMethod',obj.flags.dphiMethod);
            
            %3\ Process the image
            nSrc    = numel(obj.trs.src);
            obj.psf = psfStats.empty(nSrc,0);
            for iSrc = 1:nSrc
                obj.psf(iSrc)  = psfStats(obj.rec_(:,:,iSrc),obj.trs.tel.pupil,obj.trs.cam.wavelength,obj.trs.cam.samp,...
                    obj.trs.cam.pixelScale,obj.trs.src(iSrc),'im_ref',obj.trs.cam.image(:,:,iSrc),'fitbg',obj.flags.fitbg,'psfResolution',size(obj.trs.cam.image,1));
            end
            %4\ Prompt messages
            fprintf('-----------------------------------------\n');
            fprintf('Time for reconstructing the PSF :\t %.3g s\n', toc(tpsfr));
            fprintf('-----------------------------------------\n');
        end
    end
end




