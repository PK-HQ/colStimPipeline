function addSkippedTicks(minLim,maxLim,intervalSize,axisToAdd)
switch axisToAdd
    case {'x'}
        xlim([minLim maxLim])
        xticks(minLim:intervalSize:maxLim)
        ax=gca;
        labels = string(ax.XAxis.TickLabels); % extract
        labels(2:2:end) = nan; % remove every other one
        ax.XAxis.TickLabels = labels; % set
        xlim([minLim maxLim])
        xtickangle(0);
    case {'y'}
        ylim([minLim maxLim])
        yticks(minLim:intervalSize:maxLim)
        ax=gca;
        labels = string(ax.YAxis.TickLabels); % extract
        labels(2:2:end) = nan; % remove every other one
        ax.YAxis.TickLabels = labels; % set
        ylim([minLim maxLim])
        ytickangle(0)
end
        
end