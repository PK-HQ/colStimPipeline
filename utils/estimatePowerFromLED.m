function estimatedpower=estimatePowerFromLED(currentBlockStruct, plotflag)% given data
% function to estimate power measurement based on led percent
% Convert the string dates to datetime objects for comparison
currentDate = datetime(currentBlockStruct.date, 'InputFormat', 'yyyyMMdd');
referenceDate = datetime('20230209', 'InputFormat', 'yyyyMMdd');

% Check if currentDate is after referenceDate
if currentDate > referenceDate
    % Perform action X if currentDate is after referenceDate
    chamberid='L';
    led_percent = [5 10 15 20 25 30 40 50 60 70 75 80 90 100];
    powers = [25.3 32.7 40.2 47.7 54.8 62.6 77.6 92.2 106.0 120.3 127.0 133.5 146.3 159.0];
    area = 1080 * 1920 * .0054^2; % x y pixSize
    powerDensities = powers ./ area;
    ledpercenttoestimate=30; % 30OD16?
else
    % Perform action Y otherwise
    chamberid='R';
    led_percent = [5 10 15 20 25 30 40 50 60 70 75 80 90 100];
    powers = [25.3 32.7 40.2 47.7 54.8 62.6 77.6 92.2 106.0 120.3 127.0 133.5 146.3 159.0];
    area = 1080 * 1920 * .0054^2; % x y pixSize
    powerDensities = powers ./ area;
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
    y1 = powerDensities(1:breakpointindex);
    
    % segment 2 data
    x2 = led_percent(breakpointindex+1:end);
    y2 = powerDensities(breakpointindex+1:end);
    
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
    plot(led_percent, powerDensities, 'ko', 'markerfacecolor', 'k'); % original data points
    
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
