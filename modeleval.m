%% Evaluate neuropsych model parameters
% General params
saveflag=0;

% Define a struct array with multiple paramsStructs
paramsArray = struct(...
    'params1', struct('alphaCon', nlinspace(0, 1, 10, 'nonlinear'), 'alphaIncon', 0.001, 'optoV', 0, 'nPower', 3), ...
    'params2', struct( 'alphaCon', .01, 'betaCon', 1, 'betaIncon', .0001, 'optoV', nlinspace(0, 1, 10, 'nonlinear'), 'optoH', 0, 'optoBL',  0.00001, 'nPower', 3),...
    'params3', struct( 'alphaCon', .01, 'betaCon',1, 'betaIncon',0, 'optoV', .3, 'optoH', 0, 'optoBL', nlinspace(0, .001, 10, 'nonlinear'), 'nPower', 3));

% Determine how many sets of parameters there are
numParams = numel(fieldnames(paramsArray));

% Loop through each paramsStruct in paramsArray
for i = numParams
    % Extract the current paramsStruct by dynamically accessing fields of paramsArray
    currentParamsStruct = paramsArray.(sprintf('params%d', i));
    
    % Setup figure and subplot axes for the current paramsStruct
    figureName = sprintf('Parameter eval - Set %d', i);
    figure('Name', figureName)
    nRows = 1; nCols = 2;
    gap = .05; marginV = .01; marginH = .01;
    [hAx, ~] = tight_subplot(nRows, nCols, [gap gap], [marginV+.05 marginV+.2], [marginH marginH]);
    
    % Call the simulation function with the current paramsStruct
    simulateParameterEffects(currentParamsStruct, hAx, 'full');
    
    % Set common labels and titles, if applicable
    upFontSize(24, 0.005);

    switch saveflag
        case 1
            export_fig('Y:/users/PK/Eyal/meetings/psychmodel/modeleval.pdf','-pdf','-nocrop');
            copyfile('Y:/users/PK/Eyal/meetings/psychmodel/modeleval.pdf', 'Y:/Chip/Meta/psychmodel/modeleval.pdf')
    end
end

