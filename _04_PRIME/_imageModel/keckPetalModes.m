classdef keckPetalModes < handle
    
    properties
        resolution;
        modes;
        nModes
        pupil;
        pistonGroups;
        % Keck pupil segmented
        pupSeg;
    end        
    
    methods
        
        function obj = keckPetalModes(varargin) % constructor
            inputs = inputParser;
            inputs.addParameter('resolution',200,@isnumeric);
            inputs.addParameter('spiders','off',@ischar);
            inputs.addParameter('cobs','off',@ischar);
            inputs.addParameter('petal',true,@islogical);
            inputs.parse(varargin{:});
            obj.resolution = inputs.Results.resolution;
            
            % Construct the Keck pupil
            obj.pupSeg = obj.constructKeckPupil(inputs.Results.spiders,inputs.Results.cobs);
            obj.pupil = abs(obj.pupSeg.matrix);
            
            if inputs.Results.petal
                % Define the petal modes
                obj.pistonGroups    = obj.keckSegmentGroups();
                obj.nModes          = numel(obj.pistonGroups);
                obj.modes           = obj.keckCreateModes();
            else
                %define the zonalModes
                obj.nModes = 36;
                obj.pistonGroups = {1,2,3,4,5,6,7,8,9,...
                    10,11,12,13,14,15,16,17,18,19,...
                    20,21,22,23,24,25,26,27,28,29,...
                    30,31,32,33,34,35,36};
                obj.modes = obj.keckCreateModes();
            end
            
        end
        
        function modes = keckCreateModes(obj)
           
            modes = zeros(obj.resolution^2,obj.nModes);
            
            for k=1:obj.nModes
                % Reinit the phase
                obj.pupSeg.applyPhaseError(1:obj.pupSeg.nSegments,zeros(obj.pupSeg.nSegments,1));
                
                % Creating modes
                group = obj.pistonGroups{k};
                if ~iscell(group)
                    group = {group};
                end
                nG = numel(group);
                coef = linspace(-1,1,nG);
                for j=1:nG
                    nP = numel(group{j});
                    obj.pupSeg.applyPhaseError(group{j},coef(j)*ones(nP,1));
                end
                
                % Concatenating modes
                tmp = imag(obj.pupSeg.resize(obj.resolution));
                modes(:,k) = tmp(:)/max(tmp(:));
                
            end
        end
        
    end
    
    
    methods (Static)
        
        
        function pupSeg =  constructKeckPupil(spiders,cobs)            
            % Telescope definition
            D = 10;
            if strcmp(cobs,'on')
                cobs = 0.285;
            else
                cobs = 0;
            end

            % Segment definition
            segNSides = 6;
            segRadius = 0.9;
            segNPixels = 200;
            seg = segment(segNSides,segRadius,segNPixels);
            SVpx=load('Keck_SVpx_200.txt');

            % Pupil assembly   
            pupSeg = pupilGenerator('segRef',seg,'segCoord',SVpx,...
                'D',D,'obstructionRatio',cobs ,'unit','px','flagNoGap',1);
            
            %Spiders definition
            if strcmp(spiders,'on')
                pupilSpider.n  = 3;
                pupilSpider.angle  = [0 pi/3 2*pi/3];
                pupilSpider.width       = 0.1;
                pupilSpider.usrDefined  = 0;
                pupilSpider.usrMask     = [];
                pupSeg.spiders = pupilSpider;
                pupSeg.mkSpiders(pupilSpider);
            end
        end
        
        function groups = keckSegmentGroups()
            
            groups = cell(1,6);

            %%%%%%%%%%%%%%%%%%%%%%%%
            %          22          %
            %       23    21       %
            %    24    09    20    %
            % 25    10    08    19 %
            %    11    02    07    %
            % 26    03    01    36 %
            %    12          18    %
            % 27    04    06    35 %
            %    13    05    17    %
            % 28    14    16    34 %
            %    29    15    33    %
            %       30    32       %
            %          31          %
            %%%%%%%%%%%%%%%%%%%%%%%%
          
            %1 Stair mode y-axis
            c1 = [25 26 27 28 24 11 12 13 29 23 10 3 4 14 30];%-
            c2 = [22 9 2 5 15 31]; %0
            c3 = [21 8 1 6 16 32 20 7 18 17 33 19 36 35 34];%+
            groups{1} = {c1,c2,c3};
            
            %2 Stair mode y-axis
            c1 = [22 23 21 24 9 20 25 10 8 19 11 2 7 26 3 1 36];%-
            c2 = [12 18]; %0
            c3 = [27 4 6 35 13 5 17 28 14 16 34 29 15 33 30 32 31];%+
            groups{2} = {c1,c2,c3};
                       
            %3 Seg astig
            c1 = [22 23 21 9 2 28 27 13 29 4 6 17 33 34 35];%-
            c2 = [14 16 12 18 10 8]; %0
            c3 = [25 26 24 11 3 1 7 20 19 36 5 15 30 32 31];%+
            groups{3} = {c1,c2,c3};
            
            %4 Islands
            c1 = [25 26 24 11 3];% --
            c2 = [12 10]; %-
            c3 = [1 2 4 5 7 8 9 13 14 15 29 20 21 22 23 27 28 29 30 31 32 36]; % 0
            c4 = [16 18]; % +
            c5 = [34 35 17 33 6]; % ++
            groups{4} = {c1,c2,c3,c4,c5};
                        
            %5 Staircase modes 45 deg
            c1=[22 23 24 25];
            c2=[21 9 10 11 26];
            c3=[20 8 2 3 12 27];
            c4=[19 7 1 4 13 28];
            c5=[36 18 6 5 14 29];
            c6=[35 17 16 15 30];
            c7=[34 33 32 31];
            groups{5} = {c1,c2,c3,c4,c5,c6,c7};

            %c1=[25 26 27 28];
            %c2=[24 11 12 13 29];
            %c3=[23 10 3 4 14 30];
            %c4=[22 9 2 5 15 31];
            %c5=[21 8  1 6 16 32];
            %c6=[20 7 18 17 33];
            %c7=[19 36 35 34];
            %groups{5} = {c1,c2,c3,c4,c5,c6,c7};
                    
            
            %6 Longitudinal staircase mode
            c1=[22];
            c2=[23 21];
            c3=[24 09 20];
            c4=[25 10 8 19];
            c5=[11 2 7];
            c6=[26 3 1 36];
            c7=[12 18];
            c8=[27 4 6 35];
            c9=[13 5 17];
            c10=[28 14 16 34];
            c11=[29 15 33];
            c12=[30 32];
            c13 = [31];
            groups{6} = {c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13};
                                             
        end
        
    end
end

