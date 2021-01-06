function out = removeHotPixelsMap(frame,varargin)
inputs = inputParser;
inputs.addRequired('frame', @isnumeric);
inputs.addParameter('segmentSize', 128,@isnumeric);
inputs.addParameter('threshold', 4,@isnumeric);
inputs.parse(frame,varargin{:});

%1\ Parsing inputs
segmentSize = inputs.Results.segmentSize;
threshold = inputs.Results.threshold;
nFrame = size(frame);

%3\ Loop on segment
out = 0*frame;
nSeg = nFrame/segmentSize;
g = 0;
epsmx = 1;
epsmn = 1;

for k=1:nSeg(1)
    id_k = (k-1)*segmentSize+1:k*segmentSize;
    for j=1:nSeg(2)
        id_j = (j-1)*segmentSize+1:j*segmentSize;
        %select the segment
        seg = frame(id_k,id_j);
        %get min and max value
        while epsmx && epsmn && g<100
            [mx_seg,mn_seg] = puakoTools.minmax(seg);
            pos_mx = find(seg == mx_seg);
            pos_mn = find(seg == mn_seg);
            % sort the segment values
            seg_sort = sort(seg);
            % Thresholding
            mseg    = mean(seg_sort(2:numel(seg_sort)/2-1));
            stdseg  = std(seg_sort(2:numel(seg_sort)/2-1));
            if mx_seg - mseg >= threshold*stdseg
                seg(pos_mx(1)) = mseg;
            else
                epsmx = 0;
            end
            if mseg - mn_seg >= threshold*stdseg
                seg(pos_mn(1)) = mseg;
            else
                epsmn = 0;
            end
            g = g+1;
        end
        out(id_k,id_j) = seg;
    end
end
