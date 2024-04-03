% Assume behavioralData is your struct containing the 'fit' field

% Define the line numbers
lineNumbers = 1:3;
lineStr={'Baseline', 'Con-opto', 'Incon-opto'};
lineColors={[0.2 0.2 0.2],[0.9294, 0.1098, 0.1373]*1.05,[0, 0.0941, 0.6627]*1.25,[0.4471, 0.0353, 0.7176]*1.1};


% Iterate through each field
fields = {'rmax', 'c50', 'beta', 'exponent'};
for fieldIndex = 1:length(fields)
    % Create a new figure for the current field
    figure('Name', fields{fieldIndex});
    currentValuesLines=[];
    % Iterate through line numbers
    for lineIndex = 1:length(lineNumbers)
        % Extract values for the current field and line number
        currentValues = vertcat(behavioralData.fit(lineNumbers(lineIndex), 3, :).(fields{fieldIndex}));
        
        % Set binwidth, xlims
        if max(currentValues) < 10
            binWidth=1;
            xlimits=[0 10];
        elseif max(currentValues) >= 10 & max(currentValues) < 50
            binWidth=5;
            xlimits=[0 50];
        elseif max(currentValues) > 50
            binWidth=10;
            xlimits=[0 100];
        end
        
        % Create a subplot for the current line number
        subplot(1, 3, lineIndex); 
        
        % Plot histogram
        histogram(currentValues,'Normalization','probability','BinWidth',binWidth/2,'LineWidth',2,...
            'FaceColor',lineColors{lineIndex},'FaceAlpha',.8);
        
        % Add title for the current subplot
        % Calculate the median and MAD
        medianValue = median(currentValues);
        madValue = mad(currentValues, 1); % Use 1 for normalization
        
        behavioralData.initParams(fieldIndex,lineNumbers(lineIndex))=medianValue;
        
        % Append spread information to the title string
        spreadInfo = sprintf('\\pm %s', num2str(madValue, 1));
        titleString = sprintf('%s\nMedian \\pm MAD=%s, %s', lineStr{lineIndex}, num2str(medianValue, 3), spreadInfo);
        
        % Set the title
        title(titleString);
        
        currentValuesLines=[currentValuesLines; currentValues];
        
        xlim(xlimits)
        ylim([0 1])
        upFontSize()
        axis square

    end

end
