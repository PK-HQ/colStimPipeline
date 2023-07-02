function addSkippedTicks(minLim,maxLim,intervalSize,axisToAdd)
switch axisToAdd
    case {'x'}
        xticks(minLim:intervalSize:maxLim)
        ax=gca;
        labels = string(ax.XAxis.TickLabels); % extract
        labels(2:2:end) = nan; % remove every other one
        ax.XAxis.TickLabels = labels; % set
    case {'y'}
        yticks(minLim:intervalSize:maxLim)
        ax=gca;
        labels = string(ax.YAxis.TickLabels); % extract
        labels(2:2:end) = nan; % remove every other one
        ax.YAxis.TickLabels = labels; % set
end
        
end