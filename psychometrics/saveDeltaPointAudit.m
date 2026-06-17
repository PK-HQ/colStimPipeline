function paths = saveDeltaPointAudit(audit, outputBase)
% Save side-delta contrast-pairing audit as MAT and CSV.

    if ~istable(audit)
        error('saveDeltaPointAudit:InvalidAudit', ...
            'audit must be a table.');
    end
    outputBase = char(outputBase);
    outputFolder = fileparts(outputBase);
    if ~isempty(outputFolder) && ~isfolder(outputFolder)
        mkdir(outputFolder);
    end

    matPath = [outputBase '.mat'];
    csvPath = [outputBase '.csv'];
    deltaPointAudit = audit;
    save(matPath, 'deltaPointAudit');
    writetable(audit, csvPath);
    printDeltaAuditSaveSummary(audit);
    paths = struct('mat', matPath, 'csv', csvPath);
    fprintf('Delta point audit saved:\n  %s\n  %s\n', matPath, csvPath);
end

function printDeltaAuditSaveSummary(audit)
    if isempty(audit) || height(audit) == 0
        warning('saveDeltaPointAudit:EmptyAudit', ...
            'Delta point audit is empty; no side-specific delta pairs were saved.');
        return;
    end

    experimentIDs = unique(audit.experimentID(strlength(audit.experimentID) > 0));
    plottedRows = audit.isPlotted;
    skippedRows = ~audit.isPlotted;
    exactRows = plottedRows & contains(audit.pairingMethod, 'exact');
    toleranceRows = plottedRows & contains(audit.pairingMethod, 'tolerance');
    nearestRows = plottedRows & contains(audit.pairingMethod, 'nearest');
    nearRows = skippedRows & audit.nearMatchOutsideDefaultTolerance;
    missingZeroRows = skippedRows & contains(audit.reasonExcluded, ...
        'zero contrast unavailable');
    unmatchedNonzeroRows = skippedRows & contains(audit.reasonExcluded, ...
        'no one-to-one match within');

    fprintf(['Delta point audit summary: %d experiments | %d plotted rows | ' ...
        '%d skipped rows | exact=%d | tolerance/rank=%d | ' ...
        'nearest/tolerance=%d | missing zero=%d | unmatched nonzero=%d | ' ...
        'near match outside default tolerance=%d\n'], ...
        numel(experimentIDs), sum(plottedRows), sum(skippedRows), ...
        sum(exactRows), sum(toleranceRows), sum(nearestRows), ...
        sum(missingZeroRows), sum(unmatchedNonzeroRows), sum(nearRows));
    printDeltaAuditByPanel(audit);
    printSparseDeltaPanels(audit);
end

function printDeltaAuditByPanel(audit)
    sides = ["Horizontal", "Vertical"];
    metrics = ["deltaBias", "deltaMask"];
    for sideIdx = 1:numel(sides)
        for metricIdx = 1:numel(metrics)
            rows = audit.visualSide == sides(sideIdx) & ...
                audit.deltaMetric == metrics(metricIdx);
            if ~any(rows)
                continue;
            end
            plotted = rows & audit.isPlotted;
            exactRows = plotted & contains(audit.pairingMethod, 'exact');
            toleranceRows = plotted & contains(audit.pairingMethod, ...
                'tolerance');
            missingZeroRows = rows & ~audit.isPlotted & ...
                contains(audit.reasonExcluded, 'zero contrast unavailable');
            unmatchedRows = rows & ~audit.isPlotted & ...
                contains(audit.reasonExcluded, ...
                'no one-to-one match within');
            nearRows = rows & ~audit.isPlotted & ...
                audit.nearMatchOutsideDefaultTolerance;
            fprintf(['  %s %s: expected=%d plotted=%d exact=%d ' ...
                'tolerance=%d missingZero=%d unmatchedNonzero=%d ' ...
                'near5to10=%d\n'], ...
                sides(sideIdx), metrics(metricIdx), sum(plotted), ...
                sum(plotted), sum(exactRows), sum(toleranceRows), ...
                sum(missingZeroRows), sum(unmatchedRows), sum(nearRows));
        end
    end
end

function printSparseDeltaPanels(audit)
    [groupID, experimentGroup, sideGroup, metricGroup] = findgroups( ...
        audit.experimentID, audit.visualSide, audit.deltaMetric);
    nGroups = max(groupID);
    for groupIdx = 1:nGroups
        rows = groupID == groupIdx;
        plottedCount = sum(audit.isPlotted(rows));
        if plottedCount >= 4
            continue;
        end
        reasonRows = rows & ~audit.isPlotted & ...
            strlength(audit.reasonExcluded) > 0;
        reasons = unique(audit.reasonExcluded(reasonRows), 'stable');
        if isempty(reasons)
            reasonText = "no additional skipped audit rows";
        else
            reasonText = strjoin(reasons, '; ');
        end
        fprintf(['  fewerThan4: %s %s %s plotted=%d | %s\n'], ...
            experimentGroup(groupIdx), sideGroup(groupIdx), ...
            metricGroup(groupIdx), plottedCount, reasonText);
    end
end
