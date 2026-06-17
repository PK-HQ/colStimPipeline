function annotationInfo = addSignificanceBrackets(ax, xPairs, results, opts)
% Draw fixed-lane horizontal comparison lines and paired labels.

    if nargin < 4 || isempty(opts)
        opts = struct();
    end
    opts = applyDefaults(opts);

    currentLimits = ylim(ax);
    dataBottom = getOption(opts, 'annotationDataBottom', currentLimits(1));
    dataTop = getOption(opts, 'annotationDataTop', currentLimits(2));
    dataRange = dataTop - dataBottom;
    if ~isfinite(dataRange) || dataRange <= 0
        dataRange = max(abs(dataTop), 1);
    end

    isPercentCorrect = strcmpi(opts.annotationMode, 'percentCorrect');
    if isPercentCorrect
        dataRange = 100;
        lineSpacing = opts.percentLineSpacing;
        starLabelGap = opts.percentStarLabelGap;
        nsLabelGap = opts.percentNSLabelGap;
        firstLineY = opts.percentFirstLineY;
    else
        lineSpacing = opts.lineSpacingFraction .* dataRange;
        starLabelGap = opts.starLabelGapFraction .* dataRange;
        nsLabelGap = opts.nsLabelGapFraction .* dataRange;
        firstLineY = dataTop + opts.topPadFraction .* dataRange;
    end

    nResults = numel(results);
    drawnMask = false(nResults, 1);
    skippedReasons = repmat({''}, nResults, 1);
    pairObjects = repmat(emptyPairObject(), nResults, 1);
    highestTextY = dataTop;
    hasOneSampleComparisons = any(xPairs(:, 1) == xPairs(:, 2));

    for resultIdx = 1:nResults
        xPair = xPairs(resultIdx, :);
        lane = comparisonLane(xPair, hasOneSampleComparisons);
        lineY = firstLineY + (lane - 1) .* lineSpacing;
        if strcmp(results(resultIdx).star, 'n.s.')
            labelGap = nsLabelGap;
        else
            labelGap = starLabelGap;
        end
        textY = lineY + labelGap;
        label = results(resultIdx).star;

        pairObjects(resultIdx).x1 = xPair(1);
        pairObjects(resultIdx).x2 = xPair(2);
        pairObjects(resultIdx).lane = lane;
        pairObjects(resultIdx).yLine = lineY;
        pairObjects(resultIdx).yText = textY;
        pairObjects(resultIdx).labelGap = labelGap;
        pairObjects(resultIdx).label = label;

        if ~isfinite(results(resultIdx).adjustedP)
            skippedReasons{resultIdx} = results(resultIdx).test;
            continue;
        end
        if results(resultIdx).adjustedP >= opts.alpha && ~opts.showNS
            skippedReasons{resultIdx} = ...
                'non-significant and showNS is false';
            continue;
        end

        drawnMask(resultIdx) = true;
        highestTextY = max(highestTextY, textY);
    end

    if isPercentCorrect
        finalLimits = [0 100];
    else
        finalLimits = [currentLimits(1), ...
            highestTextY + opts.finalTopPadFraction .* dataRange];
    end
    if ~isPercentCorrect && isfield(opts, 'annotationFinalYLim') && ...
            numel(opts.annotationFinalYLim) == 2
        finalLimits = opts.annotationFinalYLim;
    end
    ylim(ax, finalLimits);

    hold(ax, 'on');
    for resultIdx = find(drawnMask(:))'
        pair = pairObjects(resultIdx);
        xLine = [pair.x1 pair.x2];
        if pair.x1 == pair.x2
            xLine = pair.x1 + [-opts.oneSampleHalfWidth ...
                opts.oneSampleHalfWidth];
        end

        pair.lineHandle = line(ax, xLine, [pair.yLine pair.yLine], ...
            'Color', 'k', ...
            'LineWidth', opts.lineWidth, ...
            'Clipping', 'off', ...
            'HandleVisibility', 'off', ...
            'Tag', 'sigLine');

        isSignificant = results(resultIdx).adjustedP < opts.alpha;
        if isSignificant
            fontSize = opts.starFontSize;
            fontWeight = 'bold';
        else
            fontSize = opts.nsFontSize;
            fontWeight = 'normal';
        end
        pair.fontSize = fontSize;
        pair.textHandle = text(ax, mean(xLine), pair.yText, pair.label, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', fontSize, ...
            'FontWeight', fontWeight, ...
            'Color', 'k', ...
            'Interpreter', 'none', ...
            'Clipping', 'off', ...
            'Tag', 'sigLabel');
        pairObjects(resultIdx) = pair;
    end

    annotationInfo = struct( ...
        'shownComparisons', find(drawnMask), ...
        'drawnMask', drawnMask, ...
        'skippedReasons', {skippedReasons}, ...
        'pairs', pairObjects, ...
        'mode', opts.annotationMode, ...
        'yDataMax', dataTop, ...
        'yRangeBase', dataRange, ...
        'lineSpacing', lineSpacing, ...
        'starLabelGap', starLabelGap, ...
        'nsLabelGap', nsLabelGap, ...
        'finalYLim', finalLimits, ...
        'nDrawn', sum(drawnMask));
end

function lane = comparisonLane(xPair, hasOneSampleComparisons)
    sortedPair = sort(xPair);
    if isequal(xPair, [1 1])
        lane = 1;
    elseif isequal(xPair, [2 2])
        lane = 2;
    elseif hasOneSampleComparisons && isequal(sortedPair, [1 2])
        lane = 3;
    elseif isequal(sortedPair, [1 2])
        lane = 1;
    elseif isequal(sortedPair, [2 3])
        lane = 2;
    elseif isequal(sortedPair, [1 3])
        lane = 3;
    else
        lane = 1;
    end
end

function pair = emptyPairObject()
    pair = struct( ...
        'x1', nan, ...
        'x2', nan, ...
        'lane', nan, ...
        'yLine', nan, ...
        'yText', nan, ...
        'labelGap', nan, ...
        'label', '', ...
        'fontSize', nan, ...
        'lineHandle', gobjects(1), ...
        'textHandle', gobjects(1));
end

function opts = applyDefaults(opts)
    defaults = struct( ...
        'alpha', 0.05, ...
        'showNS', false, ...
        'lineWidth', 1.0, ...
        'starFontSize', 13, ...
        'nsFontSize', 8, ...
        'annotationMode', 'standard', ...
        'lineSpacingFraction', 0.16, ...
        'starLabelGapFraction', 0.004, ...
        'nsLabelGapFraction', 0.010, ...
        'topPadFraction', 0.08, ...
        'finalTopPadFraction', 0.05, ...
        'percentFirstLineY', 102, ...
        'percentLineSpacing', 7, ...
        'percentStarLabelGap', 0.5, ...
        'percentNSLabelGap', 1.0, ...
        'oneSampleHalfWidth', 0.14);
    names = fieldnames(defaults);
    for idx = 1:numel(names)
        if ~isfield(opts, names{idx}) || isempty(opts.(names{idx}))
            opts.(names{idx}) = defaults.(names{idx});
        end
    end
end

function value = getOption(opts, fieldName, defaultValue)
    if isfield(opts, fieldName) && ~isempty(opts.(fieldName))
        value = opts.(fieldName);
    else
        value = defaultValue;
    end
end
