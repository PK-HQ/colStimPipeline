function paths = savePowerMetricsAudit(audit, outputBase)
% Save canonical power metrics audit as MAT and CSV.

    if ~istable(audit)
        error('savePowerMetricsAudit:InvalidAudit', ...
            'audit must be a table.');
    end
    outputBase = char(outputBase);
    outputFolder = fileparts(outputBase);
    if ~isempty(outputFolder) && ~isfolder(outputFolder)
        mkdir(outputFolder);
    end

    matPath = [outputBase '.mat'];
    csvPath = [outputBase '.csv'];
    powerMetricsAudit = audit;
    save(matPath, 'powerMetricsAudit');
    writetable(audit, csvPath);
    paths = struct('mat', matPath, 'csv', csvPath);
    fprintf('Power metrics audit saved:\n  %s\n  %s\n', matPath, csvPath);
end
