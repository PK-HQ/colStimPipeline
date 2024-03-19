function estimatedpower=estimatePowerFromLED(currentBlockStruct, plotflag)% given data
% function to estimate power measurement based on led percent
% Convert the string dates to datetime objects for comparison
currentDate = datetime(currentBlockStruct.date, 'InputFormat', 'yyyyMMdd');
referenceDate = datetime('20230209', 'InputFormat', 'yyyyMMdd');

% Check if currentDate is after referenceDate
if currentDate > referenceDate
    % Perform action X if currentDate is after referenceDate
    chamberid='L';
    led_percent = [5 10 20 30 40 50 60 70 80 90 100];
    power_measurements = [0.14 0.16 0.2 0.24 0.26 0.66 0.73 0.79 0.85 0.91 0.96];
    ledpercenttoestimate=30; % 30OD16?
else
    % Perform action Y otherwise
    chamberid='R';
    led_percent = [5 10 20 30 40 50 60 70 80 90 100];
    power_measurements = [0.14 0.16 0.2 0.24 0.26 0.66 0.73 0.79 0.85 0.91 0.96];
    ledpercenttoestimate=100; % 100OD10?
end

% initialize variables for storing the best fit
bestfit = [];
minerror = inf;
bestbreakpoint = nan;

% try breakpoints from the second to the second-last data point
for breakpointindex = 2:length(led_percent)-1
    % segment 1 data
    x1 = led_percent(1:breakpointindex);
    y1 = power_measurements(1:breakpointindex);
    
    % segment 2 data
    x2 = led_percent(breakpointindex+1:end);
    y2 = power_measurements(breakpointindex+1:end);
    
    % linear fit for segment 1
    p1 = polyfit(x1, y1, 1);
    
    % linear fit for segment 2
    p2 = polyfit(x2, y2, 1);
    
    % calculate total error
    fiterror = sum((polyval(p1, x1) - y1).^2) + sum((polyval(p2, x2) - y2).^2);
    
    % update best fit if current error is lower
    if fiterror < minerror
        minerror = fiterror;
        bestfit = {p1, p2};
        bestbreakpoint = breakpointindex;
    end
end

% plot results if a best fit was found and plotflag==1
if ~isempty(bestfit) && plotflag==1
    figure; hold on;
    plot(led_percent, power_measurements, 'ko', 'markerfacecolor', 'k'); % original data points
    
    % plot segment 1 best fit
    x1fit = linspace(min(led_percent(1:bestbreakpoint)), max(led_percent(1:bestbreakpoint)), 100);
    y1fit = polyval(bestfit{1}, x1fit);
    plot(x1fit, y1fit, 'b-', 'linewidth', 2);
    
    % plot segment 2 best fit
    x2fit = linspace(min(led_percent(bestbreakpoint+1:end)), max(led_percent(bestbreakpoint+1:end)), 100);
    y2fit = polyval(bestfit{2}, x2fit);
    plot(x2fit, y2fit, 'b-', 'linewidth', 2);
    
    xlabel('led %');
    ylabel('power measurement (mw)');
    title('two-segment linear fit to led % vs. power measurement');
    legend('data', 'linear fit', 'location', 'best');
    hold off;
    
    % display best breakpoint
    fprintf('best breakpoint at led %%: %d\n', led_percent(bestbreakpoint));
end

% estimate power for a specific led % depending on whether it falls on the
% left or right of breakpoint
if ledpercenttoestimate <= led_percent(bestbreakpoint)
    % use first segment's fit
    estimatedpower = polyval(bestfit{1}, ledpercenttoestimate);
else
    % use second segment's fit
    estimatedpower = polyval(bestfit{2}, ledpercenttoestimate);
end
end
