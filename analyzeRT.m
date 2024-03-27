function analyzeRT(TSbiasing, TSbaseline)
    % Initialize empty arrays to store time differences for consistent and inconsistent conditions
    rtCon = [];
    rtIncon = [];

    % Loop through each entry/row in the struct array
    for n = 1:numel(TSbiasing.Trial)
        % Check if fnBLK field is non-empty
        if TSbiasing.Trial(n).FlagOIBLK>0
            % Check if Graphics field length is 7 and Filename field exists
            if length(TSbiasing.Trial(n).Graphics) == 7 && isfield(TSbiasing.Trial(n).Graphics{1,3}, 'Filename')
                % Check if both filenames contain 'O00000' or 'O09000' using regex
                filename_3 = TSbiasing.Trial(n).Graphics{1,3}.Filename;
                filename_4 = TSbiasing.Trial(n).Graphics{1,4}.Filename;
                % Extract TimeSaccadeStart and TimeStimEnd values and calculate the difference
                time_diff = TSbiasing.Trial(n).TimeSaccadeStart - TSbiasing.Trial(n).TimeStimEnd;
                if contains(filename_3, 'O00000') == contains(filename_4, 'O00000') || contains(filename_3, 'O09000') == contains(filename_4, 'O09000')
                    % Store in rtCon if the time difference is positive
                    rtCon = [rtCon, time_diff];
                else
                    % Store in rtIncon if the time difference is non-positive
                    rtIncon = [rtIncon, time_diff];
                end
            end
        end
    end

    % Loop through each entry/row in the struct array
    for n = 1:numel(TSbaseline.Trial)
        % Check if fnBLK field is non-empty
        if TSbaseline.Trial(n).FlagOIBLK>0
            % Check if Graphics field length is 7 and Filename field exists
            if length(TSbaseline.Trial(n).Graphics) == 7 && isfield(TSbaseline.Trial(n).Graphics{1,3}, 'Filename')
                % Check if both filenames contain 'O00000' or 'O09000' using regex
                filename_3 = TSbaseline.Trial(n).Graphics{1,3}.Filename;
                filename_4 = TSbaseline.Trial(n).Graphics{1,4}.Filename;
                % Extract TimeSaccadeStart and TimeStimEnd values and calculate the difference
                time_diff = TSbaseline.Trial(n).TimeSaccadeStart - TSbaseline.Trial(n).TimeStimEnd;
                if contains(filename_3, 'L000') == contains(filename_4, 'O00000')
                    % Store in rtCon if the time difference is positive
                    rtH = [rtCon, time_diff];
                elseif  contains(filename_3, 'L000') == contains(filename_4, 'O09000')
                    % Store in rtIncon if the time difference is non-positive
                    rtV= [rtIncon, time_diff];
                end
            end
        end
    end

    % Plot histogram of the differences for consistent and inconsistent conditions
    plotHistogram(rtCon, rtIncon);

    % Perform statistical test and compare distributions
    compareDistributions(rtV, rtH)
    compareDistributions(rtCon, rtIncon);
end

% Plot histogram for both consistent and inconsistent conditions
function plotHistogram(rtCon, rtIncon)
    % Plot histogram of the differences for consistent conditions
    figure;
    hold on; % Hold the current plot
    histogram(rtCon,'Normalization','percentage','FaceColor','r'); % Crimson Red

    % Plot histogram of the differences for inconsistent conditions on the same axes
    histogram(rtIncon,'Normalization','percentage','FaceColor','b'); % Navy Blue
    legend('Congruent','Incongruent');
    xlabel('Time Difference (s)');
    ylabel('Frequency');
    title('Histogram of Time Differences');
    xlim([0 600])
    ylim([0 50])
    upFontSize;
end

% Compare both distributions
function compareDistributions(rtCon, rtIncon)
    % Perform the two-sample Kolmogorov-Smirnov test
    [p, h, stats] = signrank(rtCon, rtIncon);

    % Display the results
    if h == 1
        fprintf('Significantly different (p = %.4f).\n', p);
        disp(stats)
    else
        fprintf('No significantly different (p = %.4f).\n', p);
        disp(stats)
    end

    % Determine which distribution has the smaller mean
    if mean(rtCon) < mean(rtIncon)
        fprintf('rtCon has a smaller mean than rtIncon.\n');
    elseif mean(rtCon) > mean(rtIncon)
        fprintf('rtIncon has a smaller mean than rtCon.\n');
    else
        fprintf('The means of rtCon and rtIncon are equal.\n');
    end
end
