%% Evaluate neuropsych model parameters
% General params
saveflag=0;
modelComplexity='full';

% Define a struct array with multiple paramsStructs
paramsArray = struct('params1', struct('optoV', nlinspace(0, .8, 5, 'nonlinear'), 'hypothesis', 'optopercept'),...
    'params2', struct('optoV', nlinspace(0, .5, 5, 'nonlinear'), 'hypothesis', 'decisionbias'),...
    'params3', struct('optoV', nlinspace(0, .5, 5, 'nonlinear'), 'hypothesis', 'winner'));
% Determine how many sets of parameters there are
numParams = numel(fieldnames(paramsArray));

% Loop through each paramsStruct in paramsArray
for i = 1:2%1:numParams
    % Extract the current paramsStruct by dynamically accessing fields of paramsArray
    currentParamsStruct = paramsArray.(sprintf('params%d', i));
    
    % Setup figure and subplot axes for the current paramsStruct
    figureName = sprintf('Parameter eval - Set %d', i);
    hAx=figure('Name', figureName)
    nRows = 1; nCols = 1;
    gap = .05; marginV = .01; marginH = .1;
    [hAx, ~] = tight_subplot(nRows, nCols, [gap gap], [marginV+.05 marginV+.2], [marginH marginH]);
    
    % Call the simulation function with the current paramsStruct
    %simulateParameterEffectsMdl(currentParamsStruct, hAx)
    %simulateParameterEffects(currentParamsStruct, hAx, 'full');
    simulateParameterEffectsHypotheses(currentParamsStruct, hAx, modelComplexity)

    % Set common labels and titles, if applicable
    upFontSize(24, 0.005);

    switch saveflag
        case 1
            export_fig('Y:/users/PK/Eyal/meetings/psychmodel/modeleval.pdf','-pdf','-nocrop');
            copyfile('Y:/users/PK/Eyal/meetings/psychmodel/modeleval.pdf', 'Y:/Chip/Meta/psychmodel/modeleval.pdf')
    end
end

