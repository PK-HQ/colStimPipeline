%% compareModels.m
% This script compares the output curves from v_o_psychometric_3 and Optimized.
% It uses the same model parameters and contrast values.
%
% Model parameters: [a, b, l, w, g0, n, rmx, e, o]
%
lb = [0,   0,   0,   0,   0.1, 0.5,   1,  0.0,   0 ];
ub = [2,   2,   2,   2,   500, 4,   500, 1.0, 200 ];
p0 = [0.45, 0.56, 0.04, 0.59,  172, 1.9, 400, 0.39,  18];

params = [.4 .6 .1 .1,...
    150 2 200,...
    0.5 25]; %[0.1069, 0.1363, 0.0, 0.0, 60.7054, 3.548, 99.9961, 0.388, 20];
data.xBlock=[0     8    11    15    40;...
    0     7    10    13    40;...
    0     7    10    13    40];
data.yBlock=[55    40    80    95   100;...
    65    70    85   100   100;...
    40    45    65    95   100];

% Define contrast vector (using the same step as in v_o_psychometric_3)
cmin = 0; cmax = 100; cstp = 1;
%x = cmin:cstp:cmax;  % e.g. [0 5 10 ... 100]
x = sort([nlinspace(0, 100, 100, 'nonlinear')]);

% ----- Run v_o_psychometric_3 simulation -----
% v_o_psychometric_3 is a script that plots its curves.
options = struct('nTrials', 250000);
yPredicted = v_o_psychometricFn(x, params, options);
% This script should plot three curves: red (congruent), blue (incongruent), and black (control).
% ----- Set simulation options for Optimized -----
yPredictedOpt = columnarBayesMdl(x, params, options);

% The outputs are assumed to be in [0,1] (proportion correct); multiply by 100 for percent.
% Plot Optimized curves in green (for congruent) and magenta (for incongruent)
lineCol=[1 1 1; 0.9759    0.1153    0.1442; 0    0.1176    0.8284];
figure
hold on;
subplot(1,3,1)
for cond=1:3
    scatter(data.xBlock(cond,:),data.yBlock(cond,:),200, 'LineWidth', 2.5, 'Marker', 'o','markerfacecolor',lineCol(cond,:),'MarkerEdgeColor','k'); hold on
end
plot(x, yPredicted.pcntrl, 'black', 'LineWidth', 2); hold on
plot(x, yPredicted.pc, 'red', 'LineWidth', 2); hold on
plot(x, yPredicted.pic, 'blue', 'LineWidth', 2); hold on
xlabel('Gabor Contrast (%)');
ylabel('Correct (%)');
title('Bill'); upFontSize;axis square
ylim([0 100]); addSkippedTicks(0, 100, 12.5,'y')
addSkippedTicks(0, 100, 12.5,'x')

subplot(1,3,2)
for cond=1:3
    scatter(data.xBlock(cond,:),data.yBlock(cond,:),200, 'LineWidth', 2.5, 'Marker', 'o','markerfacecolor',lineCol(cond,:),'MarkerEdgeColor','k'); hold on
end
plot(x, yPredictedOpt.pcntrl, 'LineWidth', 2, 'Color', [.5 .5 .5]); hold on
plot(x, yPredictedOpt.pc, 'magenta', 'LineWidth', 2); hold on
plot(x, yPredictedOpt.pic, 'cyan', 'LineWidth', 2); hold on
title('Optimized');upFontSize
ylim([0 100]); axis square; addSkippedTicks(0, 100, 12.5,'y')
addSkippedTicks(0, 100, 12.5,'x')

subplot(1,3,3)
for cond=1:3
    scatter(data.xBlock(cond,:),data.yBlock(cond,:),200, 'LineWidth', 2.5, 'Marker', 'o','markerfacecolor',lineCol(cond,:),'MarkerEdgeColor','k'); hold on
end
plot(x, yPredicted.pcntrl, 'black', 'LineWidth', 2); hold on
plot(x, yPredicted.pc, 'red', 'LineWidth', 2); hold on
plot(x, yPredicted.pic, 'blue', 'LineWidth', 2); hold on
plot(x, yPredictedOpt.pcntrl, 'LineWidth', 2, 'Color', [.5 .5 .5]); hold on
plot(x, yPredictedOpt.pc, 'magenta', 'LineWidth', 2); hold on
plot(x, yPredictedOpt.pic, 'cyan', 'LineWidth', 2); hold on
title('Superimposed B+O');upFontSize
ylim([0 100]);axis square; addSkippedTicks(0, 100, 12.5,'y')
addSkippedTicks(0, 100, 12.5,'x')

% Add legend and labels to compare curves
legend({'Bill: Con (red)', 'Bill: Incon (blue)', 'Bill: Control (black)', ...
        'Optimized: Baseline (gray)','Optimized: Con (magenta)', 'Optimized: Incon (cyan)'}, 'Location', 'southeast','FontSize',12);
xlabel('Gabor Contrast (%)');
ylabel('Correct (%)');
hold off;
