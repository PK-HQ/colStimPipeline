function paths = saveDistributionStatsAudit(audit, outputBase)
% Save the combined distribution statistics audit as MAT and CSV.

    if ~istable(audit)
        error('saveDistributionStatsAudit:InvalidAudit', ...
            'audit must be a table.');
    end
    outputBase = char(outputBase);
    outputFolder = fileparts(outputBase);
    if ~isempty(outputFolder) && ~isfolder(outputFolder)
        mkdir(outputFolder);
    end

    matPath = [outputBase '.mat'];
    csvPath = [outputBase '.csv'];
    distributionStatsAudit = audit;
    save(matPath, 'distributionStatsAudit');
    writetable(audit, csvPath);
    paths = struct('mat', matPath, 'csv', csvPath);

    fprintf('Distribution statistics audit saved:\n  %s\n  %s\n', ...
        matPath, csvPath);
end
