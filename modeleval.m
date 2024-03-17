%% Evaluate neuropsych model parameters
% General params
saveflag=0;

% Define a struct array with multiple paramsStructs
paramsArray = struct(...
    'params1', struct('alphaCon', nlinspace(0, 1, 10, 'nonlinear'), 'alphaIncon', 0.001, 'optoV', 0, 'nPower', 4), ...
    'params2', struct('alphaIncon', nlinspace(0, 1, 10, 'nonlinear'), 'alphaCon',1, 'optoV', 0, 'nPower', 4), ...
    'params3', struct('betaCon', nlinspace(0, 1, 10, 'nonlinear'), 'betaIncon', 0, 'alphaCon', 0.1, 'alphaIncon', 0, 'optoV', 1, 'optoH', 0, 'nPower', 8), ...
    'params4', struct('betaIncon', nlinspace(0, 1, 10, 'nonlinear'), 'betaCon', 0.3, 'alphaCon', 0.1, 'alphaIncon', 0, 'optoV', 1, 'optoH', 0, 'nPower', 8), ...
    'params5', struct('betaCon', nlinspace(0, 1, 10, 'nonlinear'), 'betaIncon', 0.4, 'alphaCon', 0.1, 'alphaIncon', 0, 'optoV', 1, 'optoH', 1, 'nPower', 8), ...
    'params6', struct('betaIncon', nlinspace(0, 1, 10, 'nonlinear'), 'betaCon', 0, 'alphaCon', 0.1, 'alphaIncon', 0, 'optoV', 1, 'optoH', 1, 'nPower', 8), ...
    'params7', struct('betaCon', .3, 'betaIncon', 0.01, 'alphaCon', .1, 'alphaIncon', 0, 'optoV',  nlinspace(0, 1.5, 10, 'linear'), 'optoH', 1, 'nPower', 8));

% Determine how many sets of parameters there are
numParams = numel(fieldnames(paramsArray));

% Loop through each paramsStruct in paramsArray
for i = 7%3:numParams
    % Extract the current paramsStruct by dynamically accessing fields of paramsArray
    currentParamsStruct = paramsArray.(sprintf('params%d', i));
    
    % Setup figure and subplot axes for the current paramsStruct
    figureName = sprintf('Parameter eval - Set %d', i);
    figure('Name', figureName)
    nRows = 1; nCols = 2;
    gap = .05; marginV = .01; marginH = .01;
    [hAx, ~] = tight_subplot(nRows, nCols, [gap gap], [marginV+.05 marginV+.2], [marginH marginH]);
    
    % Call the simulation function with the current paramsStruct
    simulateParameterEffects(currentParamsStruct, hAx);
    
    % Set common labels and titles, if applicable
    upFontSize(24, 0.005);

    switch saveflag
        case 1
            export_fig('Y:/users/PK/Eyal/meetings/psychmodel/modeleval.pdf','-pdf','-nocrop');
            copyfile('Y:/users/PK/Eyal/meetings/psychmodel/modeleval.pdf', 'Y:/Chip/Meta/psychmodel/modeleval.pdf')
    end
end

