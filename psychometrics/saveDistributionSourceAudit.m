function paths = saveDistributionSourceAudit(audit, outputBase)
% Save the combined distribution source audit as MAT and CSV.

    if ~istable(audit)
        error('saveDistributionSourceAudit:InvalidAudit', ...
            'audit must be a table.');
    end
    outputBase = char(outputBase);
    outputFolder = fileparts(outputBase);
    if ~isempty(outputFolder) && ~isfolder(outputFolder)
        mkdir(outputFolder);
    end

    matPath = [outputBase '.mat'];
    csvPath = [outputBase '.csv'];
    distributionSourceAudit = audit;
    save(matPath, 'distributionSourceAudit');
    writetable(audit, csvPath);
    paths = struct('mat', matPath, 'csv', csvPath);

    fprintf('Distribution source audit saved:\n  %s\n  %s\n', ...
        matPath, csvPath);
end
