function res = getFixedParameters(x_fixed,varargin)
inputs = inputParser;
inputs.addRequired('x_fixed',@iscell);
inputs.addParameter('Cn2',1,@isnumeric);
inputs.addParameter('az',0,@isnumeric);
inputs.parse(x_fixed,varargin{:});

Cn2 = inputs.Results.Cn2;
az = inputs.Results.az;

%init if the routine fails
res.r053 = sum(Cn2);
res.Cn2 = Cn2;
res.gAO = 1;
res.gTT = 1;
res.gAl = 1;
res.az = 0;
res.map = true;
if isempty(x_fixed)
    return
end

% Checking routines
if mod(numel(x_fixed),2)
    fprintf('You must give as many values as fields for fixed parameters\n');
    return;
end
% Taking field and values
if ~isempty(x_fixed{1})
    field = x_fixed{1:2:end};
    value = x_fixed{2:2:end};
    % Check the type
    if ~ischar(field) || ~isnumeric(value)
        fprintf('You must give a fixed parameter vector with the format {''field1'',value1,''field2'',value2}\n');
        return;
    end
    
    % Fill in the structure
    if any(contains(upper(field),'R0'))
        res.r053 = value(contains(upper(field),'R0'));
    end
    if any(contains(upper(field),'CN2'))
        res.Cn2 = value(contains(upper(field),'CN2'));
    end
    if any(contains(upper(field),'GAO'))
        res.gAO = value(contains(upper(field),'GAO'));
    end
    if any(contains(upper(field),'GTT'))
        res.gTT = value(contains(upper(field),'GTT'));
    end
    if any(contains(upper(field),'GAL'))
        res.gAl = value(contains(upper(field),'GAL'));
    end
    if any(contains(upper(field),'AZ'))
        res.az = value(contains(upper(field),'AZ'));
    end
    if any(contains(upper(field),'STATICMAP'))
        res.map = value(contains(upper(field),'STATICMAP')); %must 0 or 1;
    end
   
end