function estimatedpower=estimatePowerFromLED(currentBlockStruct, bitmapData, plotflag)% given data
%% Function to estimate power measurement based on led percent

% Check whether the field `orangeLED` exists in the struct `bitmapData`
% *and* that it contains a non-empty value.  
% If both are true, use that value; otherwise, default to 100.
if isfield(bitmapData, 'orangeLED') && ~isempty(bitmapData.orangeLED)
    ledpercenttoestimate = bitmapData.orangeLED;
else
    ledpercenttoestimate = 100;
end


% Get monkey, chamber
monkeyID=str2double(currentBlockStruct.monkeyNo);
chamberL=strcmp(currentBlockStruct.chamber,'L');

% All measurements are with 580nm ND0
switch monkeyID
    case 28
        if chamberL
            % Chip L, 20240625
            led_percent = [5 10 15 20 25 30 40 50 60 70 75 80 90 100];
            powers = [25.3 32.7 40.2 47.7 54.8 62.6 77.6 92.2 106.0 120.3 127.0 133.5 146.3 159.0];
            area = 1080 * 1920 * .0054^2; % x y pixSize
            powerDensities = powers ./ area;
        else
            % Chip R, 20230208
            led_percent = [5 10 15 20 25 30 40 50 60 70 75 80 90 100];
            powers = [23.8 30.7 38 44.8 51.6 59 73.2 86.1 99.6 112.5 118 125 137.5 149];
            area = 1080 * 1920 * .0054^2; % x y pixSize
            powerDensities = powers ./ area;
        end
    case 32
        % Pepper R, 20251028
        led_percent = [5 10 15 20 25 30 40 50 60 70 75 80 90 100];
        powers = [26.5 34.2 41.8 49.4 57.1 65.0 80.1 95.0 110.2 124.8 131.0 138.0 152.0 165.5];
        area = 1080 * 1920 * .0054^2; % x y pixSize
        powerDensities = powers ./ area;
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
