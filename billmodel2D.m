% Example parameter set
clearvars;
params.a = .5;%nlinspace(0, .2, 5, 'nonlinear'); % Weight of iso-tuned visual input (Varying param 1), lifts curve but flatten, .1
params.b = .5;%nlinspace(0, 1, 5, 'nonlinear'); % Weight of ortho-tuned visual input (Varying param 2), .4
params.l = .1; % Weight of iso-tuned visual input norm, S-shape of curve, .4 to 1
params.w = .1; % Weight of ortho-tuned visual input norm, sharpening incon, any
params.e = .5; %nlinspace(.4, .6, 5, 'nonlinear'); %nlinspace(0.3, 1, 5, 'linear'); % Relative weight of optostim effects on excitation, S-shape of curve, .4 to .7
params.g0 = nlinspace(50, 100, 5, 'linear'); % Normalization constant
params.n = 6;   % Spiking exponent
params.rmx = nlinspace(0,100,5, 'linear'); % Max response
params.o = 50;  % Opto-stim effective contrast
params.ntrl = 100000; % Trials per contrast level
x = [-linspace(0, 100, 10) linspace(0, 100, 10)]; % Contrast levels

% Detect varying parameters by length > 1
param_fields = fieldnames(params);
varying_idx = find(structfun(@(v) length(v) > 1, params));
if length(varying_idx) ~= 2
    error('Exactly two parameters must be varying.');
end
param1_name = param_fields{varying_idx(1)};
param2_name = param_fields{varying_idx(2)};
param1_values = params.(param1_name);
param2_values = params.(param2_name);

% Prepare figure and subplot grid
nRows = length(param1_values);
nCols = length(param2_values);
figure;
[hAx, ~] = tight_subplot(nRows, nCols, [0.03 0.03], [0.01 0.09], [0.05 0.05]);

% Iterate over both parameters
for i = 1:nRows
    for j = 1:nCols
        % Update parameter values
        params.(param1_name) = param1_values(i);
        params.(param2_name) = param2_values(j);

        % Call normMdl
        [yOpto, yBaseline] = normMdl(x, params);

        % Determine subplot index and select axis
        subplot_idx = (i - 1) * nCols + j;
        axes(hAx(subplot_idx));

        % Plot results
        xline(0, '--', 'LineWidth', 2, 'Color', 0.3 * [1 1 1], 'HandleVisibility', 'off'); hold on;
        yline(50, '--', 'LineWidth', 2, 'Color', 0.3 * [1 1 1], 'HandleVisibility', 'off'); hold on;

        % Red for congruent
        inconIdx = 1:length(x) / 2;
        conIdx = length(x) / 2 + 1:length(x);
        plot(x(conIdx), yOpto(conIdx), 'r', 'LineWidth', 2); hold on;
        plot(x(inconIdx), yOpto(inconIdx), 'b', 'LineWidth', 2); hold on;

        % Black for baseline
        plot(x(conIdx), yBaseline(conIdx), 'k', 'LineWidth', 2); hold on;
        plot(x(inconIdx), yBaseline(inconIdx), 'k', 'LineWidth', 2); hold on;

        % Set axis properties
        xlim([-100 100]);
        ylim([0 100]); yticks(0:25:100);

        if i==1 & j==1
            % Format axes and labels
            ax = gca;
            xticks(-100:25:100); xtickangle(0);
            xx = -100:25:100;
            cellArray = cellstr(num2str(xx(:))); cellArray(2:2:end) = {''};
            xticklabels(cellArray);
            yyaxis left; ylim([0 100]); yticks(0:25:100); set(gca, 'ycolor', 'k');
            ax.YTickLabel = flipud(ax.YTickLabel);
            yyaxis right; ylim([0 100]); yticks(0:25:100); set(gca, 'ycolor', 'k');
            yyaxis left; ylabel('Correct, incongruent (%)');
            yyaxis right; ylabel('Correct, congruent (%)');
            xlabel('Gabor contrast (%)');
        else
            % Format axes and labels
            ax = gca;
            xticks(-100:25:100); xtickangle(0);
            xx = -100:25:100;
            cellArray = cellstr(num2str(xx(:))); cellArray(2:2:end) = {''};
            xticklabels([]);
            yyaxis left; ylim([0 100]); yticks(0:25:100); set(gca, 'ycolor', 'k');yticklabels([]);
            ax.YTickLabel = flipud(ax.YTickLabel);
            yyaxis right; ylim([0 100]); yticks(0:25:100); set(gca, 'ycolor', 'k');yticklabels([]);
        end
        
        % Format subplot
        xlim([-100 100]);
        ylim([0 100]);

        axis square;
        upFontSize(12, 0.02);

        % Add title for each subplot
        title(sprintf('%s=%.2f, %s=%.2f', param1_name, param1_values(i), param2_name, param2_values(j)), ...
            'FontWeight', 'normal');
    end
end

% Add labels to outer axes
suplabel({['Psychometric Functions for varying ' param1_name ' and ' param2_name], ['\Delta' param2_name]}, 't', [.1 .1 .8 .86]);
suplabel(['\Delta' param1_name], 'y');

% Save the figure
saveFlag = 1;
if saveFlag
    savefilename = sprintf('Y:/Chip/Meta/psychometrics/model/delta_%s_%s', param1_name, param2_name);
    print(gcf, [savefilename '.png'], '-dpng', '-r600');
end

function [yOpto, yBaseline] = normMdl(x, params)
    % Visual + Opto-stim Psychometric Functions
    % Inputs:
    %   x - Stimulus contrast levels
    %   params - Structure containing model parameters
    % Outputs:
    %   yOpto - Correct percentages for congruent and incongruent optostim
    %   yBaseline - Correct percentages for baseline
    
    %% Parameters
    % excitatory parameters
    a = params.a; % Weight of iso-tuned visual input
    b = params.b; % Weight of ortho-tuned visual input

    % normalization parameters
    l = params.l; % Weight of iso-tuned visual input norm
    w = params.w; % Weight of ortho-tuned visual input norm
    g0 = params.g0; % Normalization constant

    % other parameters
    n = params.n;   % Spiking exponent
    rmx = params.rmx; % Max response
    e = params.e; % Relative weight of optostim effects on excitation
    o = params.o;  % Opto-stim effective contrast

    % stimulation setup
    nTrials = 100000; % Trials per level of stimulus contrast
    contrasts = x; % Stimulus contrast levels
    nContrasts = numel(x);

    %% Storage for results
    normresp_cH_vH_oV = zeros(1, nContrasts); 
    normresp_cV_vH_oV = zeros(1, nContrasts); 
    normresp_cH_vH_oH = zeros(1, nContrasts); 
    normresp_cV_vH_oH = zeros(1, nContrasts);
    normresp_cH_vV_oV = zeros(1, nContrasts); 
    normresp_cV_vV_oV = zeros(1, nContrasts); 
    normresp_cH_vV_oH = zeros(1, nContrasts); 
    normresp_cV_vV_oH = zeros(1, nContrasts);
    
    % Calculate results for the current parameter set
    for i = 1:nContrasts
        % Excitatory and normalization signals
        resp_cH_vH_oV = contrasts(i) + b * (1 - e) * o;
        norm_cH_vH_oV = contrasts(i) + w * e * o;
        normresp_cH_vH_oV(i) = rmx * (resp_cH_vH_oV^n) / (norm_cH_vH_oV^n + g0^n);
    
        resp_cV_vH_oV = a * contrasts(i) + (1 - e) * o;
        norm_cV_vH_oV = l * contrasts(i) + e * o;
        normresp_cV_vH_oV(i) = rmx * (resp_cV_vH_oV^n) / (norm_cV_vH_oV^n + g0^n);
    
        resp_cH_vH_oH = contrasts(i) + (1 - e) * o;
        norm_cH_vH_oH = contrasts(i) + e * o;
        normresp_cH_vH_oH(i) = rmx * (resp_cH_vH_oH^n) / (norm_cH_vH_oH^n + g0^n);
    
        resp_cV_vH_oH = a * contrasts(i) + b * (1 - e) * o;
        norm_cV_vH_oH = l * contrasts(i) + w * e * o;
        normresp_cV_vH_oH(i) = rmx * (resp_cV_vH_oH^n) / (norm_cV_vH_oH^n + g0^n);
    
        resp_cH_vV_oV = a * contrasts(i) + b * (1 - e) * o;
        norm_cH_vV_oV = l * contrasts(i) + w * e * o;
        normresp_cH_vV_oV(i) = rmx * (resp_cH_vV_oV^n) / (norm_cH_vV_oV^n + g0^n);
    
        resp_cV_vV_oV = contrasts(i) + (1 - e) * o;
        norm_cV_vV_oV = contrasts(i) + e * o;
        normresp_cV_vV_oV(i) = rmx * (resp_cV_vV_oV^n) / (norm_cV_vV_oV^n + g0^n);
    
        resp_cH_vV_oH = a * contrasts(i) + (1 - e) * o;
        norm_cH_vV_oH = l * contrasts(i) + e * o;
        normresp_cH_vV_oH(i) = rmx * (resp_cH_vV_oH^n) / (norm_cH_vV_oH^n + g0^n);
    
        resp_cV_vV_oH = contrasts(i) + b * (1 - e) * o;
        norm_cV_vV_oH = contrasts(i) + w * e * o;
        normresp_cV_vV_oH(i) = rmx * (resp_cV_vV_oH^n) / (norm_cV_vV_oH^n + g0^n);
    end

    %% Simulate trials and compute performance
    trials = zeros(nTrials, 4);
    for trial = 1:nTrials
        % Generate random input and contrast index
        inp = randi([0 1], 1, 2, 'single'); % Binary inputs for trials
        contrastIdx = randi([1 nContrasts]); % Random contrast index
        trials(trial, 1:2) = inp;
        trials(trial, 3) = contrastIdx;
    
        % Determine response based on input condition
        if inp(1) == 0 && inp(2) == 0
            % Vertical column + vertical optostim
            respH = randn() + normresp_cH_vV_oV(contrastIdx);
            respV = randn() + normresp_cV_vV_oV(contrastIdx);
        elseif inp(1) == 0 && inp(2) == 1
            % Vertical column + horizontal optostim
            respH = randn() + normresp_cH_vV_oH(contrastIdx);
            respV = randn() + normresp_cV_vV_oH(contrastIdx);
        elseif inp(1) == 1 && inp(2) == 0
            % Horizontal column + vertical optostim
            respH = randn() + normresp_cH_vH_oV(contrastIdx);
            respV = randn() + normresp_cV_vH_oV(contrastIdx);
        else
            % Horizontal column + horizontal optostim
            respH = randn() + normresp_cH_vH_oH(contrastIdx);
            respV = randn() + normresp_cV_vH_oH(contrastIdx);
        end
    
        % Calculate likelihood ratio
        lr = sum(exp(-0.5 * ((respH - [normresp_cH_vH_oV, normresp_cH_vH_oH]).^2 + ...
                             (respV - [normresp_cV_vH_oV, normresp_cV_vH_oH]).^2))) / ...
             sum(exp(-0.5 * ((respH - [normresp_cH_vV_oV, normresp_cH_vV_oH]).^2 + ...
                             (respV - [normresp_cV_vV_oV, normresp_cV_vV_oH]).^2)));
    
        % Determine trial outcome
        trials(trial, 4) = double(lr > 1.0 || (lr == 1.0 && rand() > 0.5));
    end


    %% Initialize output as a structure
    output = struct(...
        'stimulusContrast', zeros(nContrasts, 1), ...
        'effectiveOptoContrast', zeros(nContrasts, 1), ...
        'congruentTrialCount', zeros(nContrasts, 1), ...
        'incongruentHVTrialCount', zeros(nContrasts, 1), ...
        'incongruentVHTrialCount', zeros(nContrasts, 1), ...
        'baselineTrialCount', zeros(nContrasts, 1), ...
        'correctCongruentTrialCount', zeros(nContrasts, 1), ...
        'correctIncongruentHVTrialCount', zeros(nContrasts, 1), ...
        'correctIncongruentVHTrialCount', zeros(nContrasts, 1), ...
        'correctBaselineTrialCount', zeros(nContrasts, 1));

    %% Process trials
    for trial = 1:nTrials
        contrastIdx = trials(trial, 3); % Index for the current contrast level
        output.stimulusContrast(contrastIdx) = contrasts(contrastIdx); % Store the contrast level
        output.effectiveOptoContrast(contrastIdx) = o; % Store optostim effective contrast
    
        % Increment counters based on input conditions
        if trials(trial, 1) == 1 && trials(trial, 2) == 1
            % Congruent condition
            output.congruentTrialCount(contrastIdx) = output.congruentTrialCount(contrastIdx) + 1; % Count congruent trials
            if trials(trial, 1) == trials(trial, 4)
                output.correctCongruentTrialCount(contrastIdx) = output.correctCongruentTrialCount(contrastIdx) + 1; % Count correct congruent trials
            end
        elseif trials(trial, 1) == 1 && trials(trial, 2) == 0
            % Horizontal-vertical incongruent condition
            output.incongruentHVTrialCount(contrastIdx) = output.incongruentHVTrialCount(contrastIdx) + 1; % Count incongruent trials (horizontal-vertical)
            if trials(trial, 1) == trials(trial, 4)
                output.correctIncongruentHVTrialCount(contrastIdx) = output.correctIncongruentHVTrialCount(contrastIdx) + 1; % Count correct incongruent trials
            end
        elseif trials(trial, 1) == 0 && trials(trial, 2) == 1
            % Vertical-horizontal incongruent condition
            output.incongruentVHTrialCount(contrastIdx) = output.incongruentVHTrialCount(contrastIdx) + 1; % Count incongruent trials (vertical-horizontal)
            if trials(trial, 1) == trials(trial, 4)
                output.correctIncongruentVHTrialCount(contrastIdx) = output.correctIncongruentVHTrialCount(contrastIdx) + 1; % Count correct incongruent trials
            end
        elseif trials(trial, 1) == 0 && trials(trial, 2) == 0
            % Baseline condition
            output.baselineTrialCount(contrastIdx) = output.baselineTrialCount(contrastIdx) + 1; % Count baseline trials
            if trials(trial, 1) == trials(trial, 4)
                output.correctBaselineTrialCount(contrastIdx) = output.correctBaselineTrialCount(contrastIdx) + 1; % Count correct baseline trials
            end
        end
    end

    
    %% Compute probabilities
    probCorrectCongruent = (output.correctCongruentTrialCount + output.correctBaselineTrialCount) ./ ...
                           (output.congruentTrialCount + output.baselineTrialCount); % total correct con+bl / total count
    
    probCorrectIncongruent = (output.correctIncongruentHVTrialCount + output.correctIncongruentVHTrialCount) ./ ...
                             (output.incongruentHVTrialCount + output.incongruentVHTrialCount); % total correct incon / total count
    
    probCorrectBaseline = (output.correctCongruentTrialCount + output.correctIncongruentHVTrialCount + ...
                           output.correctIncongruentVHTrialCount + output.correctBaselineTrialCount) ./ ...
                          (output.congruentTrialCount + output.incongruentHVTrialCount + ...
                           output.incongruentVHTrialCount + output.baselineTrialCount);



    inconIdx = 1:nContrasts / 2;
    conIdx = nContrasts / 2 + 1:nContrasts;

    yOpto = [100 - (probCorrectIncongruent(inconIdx) * 100), probCorrectCongruent(conIdx) * 100];
    yBaseline = [100 - (probCorrectBaseline(inconIdx) * 100), probCorrectBaseline(conIdx) * 100];
end

%{
function [yOpto, yBaseline] = normMdl(x, params)
    % Visual + Opto-stim Psychometric Functions
    % Inputs:
    %   x - Stimulus contrast levels
    %   params - Structure containing model parameters
    % Outputs:
    %   yOpto - Correct percentages for congruent and incongruent optostim
    %   yBaseline - Correct percentages for baseline
    
    %% Parameters
    % excitatory parameters
    a = params.a; % Weight of iso-tuned visual input
    b = params.b; % Weight of ortho-tuned visual input

    % normalization parameters
    l = params.l; % Weight of iso-tuned visual input norm
    w = params.w; % Weight of ortho-tuned visual input norm
    g0 = params.g0; % Normalization constant

    % other parameters
    n = params.n;   % Spiking exponent
    rmx = params.rmx; % Max response
    e = params.e; % Relative weight of optostim effects on excitation
    o = params.o;  % Opto-stim effective contrast

    % stimulation setup
    ntrl = params.ntrl; % Trials per level of stimulus contrast
    c = x; % Stimulus contrast levels
    nContrasts = numel(x);

    %% Storage for results
    uhhv = zeros(1, nContrasts); uvhv = zeros(1, nContrasts); 
    uhhh = zeros(1, nContrasts); uvhh = zeros(1, nContrasts);
    uhvv = zeros(1, nContrasts); uvvv = zeros(1, nContrasts); 
    uhvh = zeros(1, nContrasts); uvvh = zeros(1, nContrasts);

    % Calculate results for the current parameter set
    for i = 1:nContrasts
        % Excitatory and normalization signals
        ehhv = c(i) + b * (1 - e) * o;
        ghhv = c(i) + w * e * o;
        uhhv(i) = rmx * (ehhv^n) / (ghhv^n + g0^n);

        evhv = a * c(i) + (1 - e) * o;
        gvhv = l * c(i) + e * o;
        uvhv(i) = rmx * (evhv^n) / (gvhv^n + g0^n);

        ehhh = c(i) + (1 - e) * o;
        ghhh = c(i) + e * o;
        uhhh(i) = rmx * (ehhh^n) / (ghhh^n + g0^n);

        evhh = a * c(i) + b * (1 - e) * o;
        gvhh = l * c(i) + w * e * o;
        uvhh(i) = rmx * (evhh^n) / (gvhh^n + g0^n);

        ehvv = a * c(i) + b * (1 - e) * o;
        ghvv = l * c(i) + w * e * o;
        uhvv(i) = rmx * (ehvv^n) / (ghvv^n + g0^n);

        evvv = c(i) + (1 - e) * o;
        gvvv = c(i) + e * o;
        uvvv(i) = rmx * (evvv^n) / (gvvv^n + g0^n);

        ehvh = a * c(i) + (1 - e) * o;
        ghvh = l * c(i) + e * o;
        uhvh(i) = rmx * (ehvh^n) / (ghvh^n + g0^n);

        evvh = c(i) + b * (1 - e) * o;
        gvvh = c(i) + w * e * o;
        uvvh(i) = rmx * (evvh^n) / (gvvh^n + g0^n);
    end

    %% Simulate trials and compute performance
    trls = zeros(ntrl, 4);
    for j = 1:ntrl
        inp = randi([0 1], 1, 2, 'single'); % Binary inputs for trials
        ic = randi([1 nContrasts]);
        trls(j, 1:2) = inp;
        trls(j, 3) = ic;

        if inp(1) == 0 && inp(2) == 0
            rh = randn() + uhvv(ic);
            rv = randn() + uvvv(ic);
        elseif inp(1) == 0 && inp(2) == 1
            rh = randn() + uhvh(ic);
            rv = randn() + uvvh(ic);
        elseif inp(1) == 1 && inp(2) == 0
            rh = randn() + uhhv(ic);
            rv = randn() + uvhv(ic);
        else
            rh = randn() + uhhh(ic);
            rv = randn() + uvhh(ic);
        end

        lr = sum(exp(-0.5 * ((rh - [uhhv, uhhh]).^2 + (rv - [uvhv, uvhh]).^2))) / ...
             sum(exp(-0.5 * ((rh - [uhvv, uhvh]).^2 + (rv - [uvvv, uvvh]).^2)));
        trls(j, 4) = double(lr > 1.0 || (lr == 1.0 && rand() > 0.5));
    end

    %% Process trials
    output = zeros(nContrasts, 10); % Columns for trial-level data (e.g., congruent/incongruent counts and correct responses)
    for j = 1:ntrl
        ic = trls(j, 3); % Index for the current contrast level
        output(ic, 1) = c(ic); % Store the contrast level
        output(ic, 2) = o; % Store optostim effective contrast

        % Increment counters based on input conditions
        if trls(j, 1) == 1 && trls(j, 2) == 1
            % Congruent condition
            output(ic, 3) = output(ic, 3) + 1; % Count congruent trials
            if trls(j, 1) == trls(j, 4)
                output(ic, 7) = output(ic, 7) + 1; % Count correct congruent trials
            end
        elseif trls(j, 1) == 1 && trls(j, 2) == 0
            % Horizontal-vertical incongruent condition
            output(ic, 4) = output(ic, 4) + 1; % Count incongruent trials (horizontal-vertical)
            if trls(j, 1) == trls(j, 4)
                output(ic, 8) = output(ic, 8) + 1; % Count correct incongruent trials
            end
        elseif trls(j, 1) == 0 && trls(j, 2) == 1
            % Vertical-horizontal incongruent condition
            output(ic, 5) = output(ic, 5) + 1; % Count incongruent trials (vertical-horizontal)
            if trls(j, 1) == trls(j, 4)
                output(ic, 9) = output(ic, 9) + 1; % Count correct incongruent trials
            end
        elseif trls(j, 1) == 0 && trls(j, 2) == 0
            % Baseline condition
            output(ic, 6) = output(ic, 6) + 1; % Count baseline trials
            if trls(j, 1) == trls(j, 4)
                output(ic, 10) = output(ic, 10) + 1; % Count correct baseline trials
            end
        end
    end

    % Compute probabilities
    pc = (output(:, 7) + output(:, 10)) ./ (output(:, 3) + output(:, 6));
    pic = (output(:, 8) + output(:, 9)) ./ (output(:, 4) + output(:, 5));
    pcntrl = (output(:, 7) + output(:, 8) + output(:, 9) + output(:, 10)) ./ ...
             (output(:, 3) + output(:, 4) + output(:, 5) + output(:, 6));

    inconIdx = 1:nContrasts / 2;
    conIdx = nContrasts / 2 + 1:nContrasts;

    yOpto = [100 - (pic(inconIdx) * 100), pc(conIdx) * 100];
    yBaseline = [100 - (pcntrl(inconIdx) * 100), pcntrl(conIdx) * 100];
end
%}


