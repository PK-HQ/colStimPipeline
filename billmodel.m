% Example parameter set
clearvars;
params.a = 0.1; %Weight of iso-tuned visual input
params.b = 0.05; %Weight of ortho-tuned visual input
params.l = 0.0; %Weight of iso-tuned visual input norm
params.e = .35:.05:.8; %Relative weight of optostim effects on excitation
params.w = 0.0; %Weight of ortho-tuned visual input norm
params.g0 = 30;
params.n = 2;
params.rmx = 3;
params.o = 40;
params.ntrl = 100000;
x = [-linspace(0, 100, 10) linspace(0, 100, 10)]; %sort([-nlinspace(0, 50, 6, 'nonlinear') nlinspace(0, 50, 6, 'nonlinear')]);
% Equation
[yOpto, yBaseline] = normMdl(x, params);

function [yOpto, yBaseline] = normMdl(x, params)
% Visual + Opto-stim Psychometric Functions
% random contrast
%close all;
%clearvars;
%% Parameters
% excitatory parameters
a_values = params.a; % Range of values for alpha (visual excitatory cross-talk parameter)
b_values = params.b; % beta, opto-stim excitatory cross-talk parameter

% normalization parameters
l_values = params.l; % lambda, visual normalization cross-talk parameter
w_values = params.w; % omega, opto-stim normalization cross-talk parameter
g0_values = params.g0; % normalization constant

% other parameters
n_values = params.n;   % spiking exponent
rmx_values = params.rmx; % max response
e_values =params.e; %0.5; % eta, relative weight of optostim on normalization vs. excitation
o_values = params.o;  % opto-stim effective contrast

% stimulation setup
ntrl = params.ntrl; % trials per level of stimulus contrast
c = x; %(cmin:cstp:cmax);
nlev = numel(x);

% Detect which parameter is varying (has more than 1 element)
paramsCell = {a_values, b_values, l_values, w_values, g0_values, n_values, rmx_values, e_values, o_values};
param_names = {'a', 'b', 'l', 'w', 'g0', 'n', 'rmx', 'e', 'o'};
param_names_full = {'\alpha', '\beta', '\lambda', '\omega', 'Normalization constant', 'Spiking exponent', 'Max response', '\eta', 'Optostim strength'};
varying_param_idx = find(cellfun(@(x) length(x) > 1, paramsCell)); % Find the varying parameter

if isempty(varying_param_idx)
    error('No parameter is set to a range of values.');
end

% Get the varying parameter values and its name
varying_param = paramsCell{varying_param_idx};
varying_param_name = param_names{varying_param_idx};
varying_param_name_full = param_names_full{varying_param_idx};

fprintf('Varying parameter: %s\n', varying_param_name);

% Set up color saturation for plotting
colorSaturation = linspace(0.2, 1, length(varying_param));  % Color saturation for red (congruent)

% Create figure
figure;
nRows=2;nCols=5;
gap=0;marginV=.0;marginH=.05;
[hAx,~]=tight_subplot(nRows,nCols,[gap+.12 gap+.085], [marginV+.15 marginV+.15], [marginH marginH]);

% Loop over the varying parameter
for k = 1:length(varying_param)
    % Set the varying parameter based on the current iteration
    switch varying_param_name
        case 'a'
            a = varying_param(k);
            b = b_values; l = l_values; w = w_values; g0 = g0_values; n = n_values; rmx = rmx_values; e = e_values; o = o_values;
        case 'b'
            b = varying_param(k);
            a = a_values; l = l_values; w = w_values; g0 = g0_values; n = n_values; rmx = rmx_values; e = e_values; o = o_values;
        case 'l'
            l = varying_param(k);
            a = a_values; b = b_values; w = w_values; g0 = g0_values; n = n_values; rmx = rmx_values; e = e_values; o = o_values;
        case 'w'
            w = varying_param(k);
            a = a_values; b = b_values; l = l_values; g0 = g0_values; n = n_values; rmx = rmx_values; e = e_values; o = o_values;
        case 'g0'
            g0 = varying_param(k);
            a = a_values; b = b_values; l = l_values; w = w_values; n = n_values; rmx = rmx_values; e = e_values; o = o_values;
        case 'n'
            n = varying_param(k);
            a = a_values; b = b_values; l = l_values; w = w_values; g0 = g0_values; rmx = rmx_values; e = e_values; o = o_values;
        case 'rmx'
            rmx = varying_param(k);
            a = a_values; b = b_values; l = l_values; w = w_values; g0 = g0_values; n = n_values; e = e_values; o = o_values;
        case 'e'
            e = varying_param(k);
            a = a_values; b = b_values; l = l_values; w = w_values; g0 = g0_values; n = n_values; rmx = rmx_values; o = o_values;
        case 'o'
            o = varying_param(k);
            a = a_values; b = b_values; l = l_values; w = w_values; g0 = g0_values; n = n_values; rmx = rmx_values; e = e_values;
    end

    % Storage for results
    uhhv = zeros(1, nlev); uvhv = zeros(1, nlev); uhhh = zeros(1, nlev); uvhh = zeros(1, nlev);
    uhvv = zeros(1, nlev); uvvv = zeros(1, nlev); uhvh = zeros(1, nlev); uvvh = zeros(1, nlev);

    % Calculate results for current parameter set
    for i = 1:nlev
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

    % Simulate trials and compute performance
    trls = zeros(ntrl, 4);
    for j = 1:ntrl
        inp = randi([0 1], 1, 2, 'single'); %'logical'
        ic = randi([1 nlev]);
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

        num = sum(exp(-0.5 * ((rh - uhhv).^2 + (rv - uvhv).^2)) + exp(-0.5 * ((rh - uhhh).^2 + (rv - uvhh).^2)));
        den = sum(exp(-0.5 * ((rh - uhvv).^2 + (rv - uvvv).^2)) + exp(-0.5 * ((rh - uhvh).^2 + (rv - uvvh).^2)));
        lr = num / den;
        if lr > 1.0
            trls(j, 4) = 1;
        elseif lr == 1
            if rand() > 0.5
                trls(j, 4) = 1;
            end
        end
    end

    % Process trials
    output = zeros(nlev, 10);
    for j = 1:ntrl
        ic = trls(j, 3);
        output(ic, 1) = c(ic);
        output(ic, 2) = o;
        if trls(j, 1) == 1 && trls(j, 2) == 1
            output(ic, 3) = output(ic, 3) + 1;
            if trls(j, 1) == trls(j, 4)
                output(ic, 7) = output(ic, 7) + 1;
            end
        elseif trls(j, 1) == 1 && trls(j, 2) == 0
            output(ic, 4) = output(ic, 4) + 1;
            if trls(j, 1) == trls(j, 4)
                output(ic, 8) = output(ic, 8) + 1;
            end
        elseif trls(j, 1) == 0 && trls(j, 2) == 1
            output(ic, 5) = output(ic, 5) + 1;
            if trls(j, 1) == trls(j, 4)
                output(ic, 9) = output(ic, 9) + 1;
            end
        elseif trls(j, 1) == 0 && trls(j, 2) == 0
            output(ic, 6) = output(ic, 6) + 1;
            if trls(j, 1) == trls(j, 4)
                output(ic, 10) = output(ic, 10) + 1;
            end
        end
    end

    % Calculate performance
    pc = (output(:, 7) + output(:, 10)) ./ (output(:, 3) + output(:, 6));
    pic = (output(:, 8) + output(:, 9)) ./ (output(:, 4) + output(:, 5));
    pcntrl = (output(:,7)+output(:,8)+output(:,9)+output(:,10))...
        ./(output(:,3)+output(:,4)+output(:,5)+output(:,6));
    
    % Get pc, pic, pcntrl
    inconIdx=1:nlev/2;
    conIdx=nlev/2 + 1 : nlev;
    yOpto=[100-(pic(inconIdx) * 100) pc(conIdx) * 100];
    yBaseline=[100-(pcntrl(inconIdx)*100) pcntrl(conIdx)*100];
    
    % Plot with varying saturation
    axes(hAx(k));
    xline(0,'--','LineWidth',2,'Color',.3*[1 1 1],'HandleVisibility','off'); hold on
    yline(50,'--','LineWidth',2,'Color',.3*[1 1 1],'HandleVisibility','off'); hold on
    % Red for congruent
    plot(output(conIdx, 1), yOpto(conIdx), 'Color', [1 0 0 colorSaturation(k)], 'LineWidth', 2); hold on; 
    % Blue for incongruent
    plot(output(inconIdx, 1), yOpto(inconIdx), 'Color', [0 0 1 colorSaturation(k)], 'LineWidth', 2); hold on; % Blue for incongruent   
    % Black for baseline
    
    plot(output(conIdx,1),yBaseline(conIdx),'Color', [0 0 0 colorSaturation(k)], 'LineWidth',2); hold on;
    plot(output(inconIdx,1),yBaseline(inconIdx),'Color', [0 0 0 colorSaturation(k)], 'LineWidth',2); hold on;
    
    
   % Dual axes
    ax = gca;
    xticks([-100:25:100]); xtickangle(0)
    xx = ([-100:25:100]);  
    cellArray = cellstr(num2str(xx(:))); cellArray(2:2:end) = {''}; xticklabels(cellArray)
    yyaxis left; ylim([0 100]); yticks([0:25:100]); set(gca,'ycolor','k') 
    ax.YTickLabel = flipud(ax.YTickLabel);
    yyaxis right; ylim([0 100]); yticks([0:25:100]); set(gca,'ycolor','k') 
    % Axis labels
    xlim([-100 100]);
    ylim([0 100]); yticks([0:25:100])
    upFontSize(20,.025)
    title([varying_param_name_full ': ' num2str(varying_param(k),2)],'FontWeight','normal')
    axis square

    % Labels
    if k==1
        yyaxis left; ylim([0 100]); yticks([0:25:100]); set(gca,'ycolor','k'); ylabel('Correct, incongruent (%)'); 
        yyaxis right; ylim([0 100]); yticks([0:25:100]); set(gca,'ycolor','k'); ylabel('Correct, congruent (%)');
        %legend('congruent', 'incongruent', 'Location',);
        ylim([0 100]);xlim([-100 100]);
        xlabel('Gabor contrast (%)')
    end
end
suplabel(['Psychometric Functions with varying ', varying_param_name_full],'t',[.1 .1 .78 .83])
upFontSize(14,.025)


% Add params table
paramsCell{varying_param_idx}='Varies';
T = array2table(paramsCell,'VariableNames',param_names);
% Get the table in string form.
TString = evalc('disp(T)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','FontName',FixedWidth,'FontSize',10,'EdgeColor','none',...
    'Units','Normalized','Position',[.63 .5 1 1]);

% Save
saveFlag=0;
monkeyName='sim';
savefilename = ['Y:\Chip\Meta\psychometrics\model\delta_' varying_param_name];
% Adjust figure properties
set(gcf, 'PaperPositionMode', 'auto');  % Use the screen dimensions for printing
set(gcf, 'PaperUnits', 'inches');       % Set units to inches
set(gcf, 'PaperSize', [8, 6]);          % Set the paper size (match aspect ratio)

% Optional: Resize figure window for desired size
%set(gcf, 'Position', [100, 100, 800, 600]);  % [x, y, width, height] in pixels

% Save as PNG
print(gcf, [savefilename '.png'], '-dpng', '-r600');
end