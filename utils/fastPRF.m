blockParam = { 'OD3.0 @ 5%', 'OD3.0 @ 50%', 'OD3.0 @ 100%',...
    'OD2.0 @ 5%', 'OD2.0 @ 20%', 'OD2.0 @ 30%','OD2.0 @ 50%',...
    'OD1.6 @ 5%','OD1.6 @ 10%'};
x = 0:.05:1.15;
runID = [6 7 8 9 10 11 12 5 4];
strID = [1:numel(blockParam)];
% Define grayscale colors (values between 0 = black and 1 = white)
lineColors = slanCM('viridis',numel(runID));  

figure('Name','PRF')
colorID=1;
for i = 1:numel(runID)
    load(['Y:\Pepper\Pepper20250328\run' num2str(runID(i)) '\Data_vdaqlog.mat'])
    if i==7
        lw=4;
    else
        lw=2;
    end
    plot(x, (TCCond{1, 9}.sum-TCCond{1, 9}.sum(1)) - (TCCond{1, 1}.sum-TCCond{1, 1}.sum(1)), ...
        'LineWidth', lw, ...
        'Color', lineColors(colorID,:), ...
        'DisplayName', blockParam{strID(i)});
    hold on
    colorID=colorID+1;
end
xlim([0 1.2])
ylim([0 3])
addSkippedTicks(0, 1.15, .1,'x')
upFontSize
axis square
legend
title('Full illumination ROI timecourse','FontWeight','Normal')
ylabel('\DeltaF/F')
xlabel('Time (s)')
