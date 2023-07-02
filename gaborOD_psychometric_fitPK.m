%% gaborOD_pedestal_psychometric.m
%close all; %clear all;
for sessions=1:10
    load(['Y:\' datastruct(1).monkey '\Loki20230407\M22D20230407R8' 'TS.mat']);
    
    %
end
load('Y:\Loki\Loki20230407\M22D20230407R9TS.mat');
TS2 = TS;
load('Y:\Loki\Loki20230407\M22D20230407R10TS.mat');
TS3 = TS;

[contrast1, contrast_trial_1, choice_trial_1] = get_choice_contrastOD_HV(TS1);
[contrast2, contrast_trial_2, choice_trial_2] = get_choice_contrastOD_HV(TS2);
[contrast3, contrast_trial_3, choice_trial_3] = get_choice_contrastOD_HV(TS3);

%collate blocks
contrast = [contrast1 contrast2 contrast3];
contrast_trial = [contrast_trial_1 contrast_trial_2 contrast_trial_3];
choice_trial = [choice_trial_1 choice_trial_2 choice_trial_3]; 

choice_proportion = nan(1, length(contrast));
for i=1:length(contrast)
    n = find(contrast_trial==contrast(i)); % collects all indeces (trial #) out of selected complete trials corresponding to a given contrast at both coordinates ; L or R
    choice_contrast= choice_trial(n);  % collects the choices of the abovementioned indeces (trial#)
    choice_proportion(i) = mean(choice_contrast);
end

[contrast j] = sort(contrast);
choice_proportion = choice_proportion(j); % These functions may not matter

%% Monitor behavioral data before fitting
figure; hold on;
plot(contrast, choice_proportion, 'k.','MarkerSize',15);
ylim([0 1]);
xlabel('Stimulus contrast (H-V, %)','FontSize',12);
ylabel('P(vertical choice)','FontSize',12);
title('Loki 4/7/2023/R8, R9, R10','FontSize',15);

%% Fit PK functions & generate fitted curve
[mu,sigma,beta,logLikelihood] = fitPsyContrast(contrast, choice_proportion);
param = [mu sigma beta];

max_contrast = max(abs(contrast));
x = -max_contrast:1:max_contrast;
x(x == 0) = []; %% curve fit does not include 0 so exclude

curveFit = dualNormCDF(param, x);

%% Plot raw data and fitted curve
figure; hold on;
set(gcf,'position',[200,200,800,600]);
plot(contrast, choice_proportion, 'k.','MarkerSize',15);
ylim([0 1]);
xlabel('Stimulus contrast (H-V, %)','FontSize',12);
ylabel('P(vertical choice)','FontSize',12);
title('Loki 4/7/2023/R8, R9, R10','FontSize',15);
subtitle(['midpoint = ', num2str(mu), '; scaling = ', num2str(beta)])
plot(x, curveFit,'-r', 'Linewidth', 2);









%% Ver 2
%% Fit with lsqcurvefit
sOptOptm = ...
    optimset('LargeScale','On', ...
    'TolFun',1e-8,'TolX',1e-8, ...
    'MaxFunEvals',1e4,'MaxIter',1e4);

% Adjusted Sigmoidal: [midpoint, sigma, beta]
PI = [10, 5, 0.5];
PLB = [20   1    0];
PUB = [40   15   1];

max_contrast = max(abs(contrast));
x = -max_contrast:1:max_contrast;
x(x == 0) = [];

params = lsqcurvefit('dualNormCDF',PI, x,y,PLB,PUB,sOptOptm);
curve_fit = dualNormCDF(params, x);

%% Plot
figure; hold on;
set(gcf,'position',[200,200,800,600]);
plot(contrast, choice_proportion, 'k.','MarkerSize',15);
ylim([0 1]);
xlabel('Stimulus contrast (H-V, %)','FontSize',12);
ylabel('P(vertical choice)','FontSize',12);
title('Loki 4/7/2023/R8, R9, R10','FontSize',15);
subtitle(['midpoint = ', num2str(params(1)), '; scaling = ', num2str(params(3))])
plot(x, curve_fit,'-r', 'Linewidth', 2);