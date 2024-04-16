function compareAndPlotDistributions(paramsStruct, chamberWanted, saveFlag)
    % Filenames
    filename=['fitParamsConstrained' chamberWanted];
    monkeyName='Chip';

    % Naming and colors
    fieldNames = {'rmax', 'exponent', 'c50', 'beta'};
    lineColors = {[0.2, 0.2, 0.2], [0.9294, 0.1098, 0.1373] * 1.05, [0, 0.0941, 0.6627] * 1.25};
    legendNames = {'Baseline', 'Con-opto', 'Incon-opto'}; % Legend names for param1, param2, param3

    for fieldIndex = 1:numel(fieldNames)
        fieldName = fieldNames{fieldIndex};
        % Extract the full field data
        fieldData = vertcat(paramsStruct.(fieldName));
        
        % Prepare the data subsets
        param1 = fieldData(:,1);
        param2 = fieldData(:,2);
        param3 = fieldData(:,3);

        % Ensure vectors are the same length for pairwise comparison
        minLength = min([length(param1), length(param2), length(param3)]);
        param1 = param1(1:minLength);
        param2 = param2(1:minLength);
        param3 = param3(1:minLength);

        % Statistical tests
        [pVal12, ~, ~] = signrank(param1, param2);
        [pVal13, ~, ~] = signrank(param1, param3);
        [pVal23, ~, ~] = signrank(param2, param3);

        % Plot distributions
        figure('Name', sprintf('%s Comparisons', fieldName));
        pairs = {{param1, param2}, {param1, param3}, {param2, param3}};
        pVals = [pVal12, pVal13, pVal23];
        colorPairs = {[1, 2], [1, 3], [2, 3]};
        pairLegendNames = {{legendNames{1}, legendNames{2}}, {legendNames{1}, legendNames{3}}, {legendNames{2}, legendNames{3}}};

        for i = 1:length(pairs)
            subplot(1, 3, i);
            pair = pairs{i};
            combinedParams = [pair{1}; pair{2}];
            currentMax = max(combinedParams);
            currentMin = min(combinedParams);

            if currentMax < 10
                binWidth = 1;
                xlimits = [0, currentMax+1];
            elseif currentMax >= 10 && currentMax < 50
                binWidth = 5;
                xlimits = [0, currentMax+5];
            else
                binWidth = 10;
                xlimits = [0, 100];
            end

            colorIndex1 = colorPairs{i}(1);
            colorIndex2 = colorPairs{i}(2);

            histogram(pair{1}, 'Normalization', 'probability', 'BinWidth', binWidth, 'FaceColor', lineColors{colorIndex1}, 'FaceAlpha', 0.7);
            xline(median(pair{1}),'--','LineWidth',2,'HandleVisibility','off');
            hold on;
            histogram(pair{2}, 'Normalization', 'probability', 'BinWidth', binWidth, 'FaceColor', lineColors{colorIndex2}, 'FaceAlpha', 0.7);
            xline(median(pair{2}),'--','LineWidth',2,'HandleVisibility','off');

            xlim(xlimits);
            ylim([0 1]);
            axis square;
            title(sprintf('p = %.3f', pVals(i)));
            xlabel([upper(fieldName(1)) fieldName(2:end)]);
            ylabel('Probability');
            upFontSize()
            legend(pairLegendNames{i}, 'Location', 'northwest');
        end
        suplabel(upper(fieldName),'t',[.1 .1 .835 .78]);upFontSize()

        % Save as PDF
        savePDF(filename, monkeyName, saveFlag, fieldIndex)
    end
end
