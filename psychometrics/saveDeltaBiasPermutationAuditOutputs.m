function paths = saveDeltaBiasPermutationAuditOutputs(monkeyName, chamberWanted, modelType, clusterMdl, aggregateFits)
% Save deltaBias permutation inspection audit tables for standalone reports.

    outputDir = fullfile('Y:\users\PK\colStimPipeline', 'outputs', 'deltaBiasPermutation');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    tag = sprintf('%s_%s_%s', monkeyName, chamberWanted, modelType);

    experimentContrasts = table();
    experimentSummary = table();
    clusterContrasts = table();
    clusterSummary = table();

    if isfield(clusterMdl, 'deltaBiasPermutationExperimentContrasts')
        experimentContrasts = clusterMdl.deltaBiasPermutationExperimentContrasts;
    end
    if isfield(clusterMdl, 'deltaBiasPermutationExperimentSummary')
        experimentSummary = clusterMdl.deltaBiasPermutationExperimentSummary;
    end
    for idx = 1:numel(aggregateFits)
        if isfield(aggregateFits(idx), 'merged') && ...
                isfield(aggregateFits(idx).merged, 'deltaBiasPermutationClusterContrasts')
            clusterContrasts = [clusterContrasts; ...
                aggregateFits(idx).merged.deltaBiasPermutationClusterContrasts]; %#ok<AGROW>
        end
        if isfield(aggregateFits(idx), 'merged') && ...
                isfield(aggregateFits(idx).merged, 'deltaBiasPermutationClusterSummary')
            clusterSummary = [clusterSummary; ...
                aggregateFits(idx).merged.deltaBiasPermutationClusterSummary]; %#ok<AGROW>
        end
    end

    paths = struct();
    paths.mat = fullfile(outputDir, sprintf('deltaBiasPermutationAudit_%s.mat', tag));
    paths.experimentContrasts = fullfile(outputDir, sprintf('deltaBiasPermutationExperimentContrasts_%s.csv', tag));
    paths.experimentSummary = fullfile(outputDir, sprintf('deltaBiasPermutationExperimentSummary_%s.csv', tag));
    paths.clusterContrasts = fullfile(outputDir, sprintf('deltaBiasPermutationClusterContrasts_%s.csv', tag));
    paths.clusterSummary = fullfile(outputDir, sprintf('deltaBiasPermutationClusterSummary_%s.csv', tag));

    save(paths.mat, 'experimentContrasts', 'experimentSummary', ...
        'clusterContrasts', 'clusterSummary');
    writetable(experimentContrasts, paths.experimentContrasts);
    writetable(experimentSummary, paths.experimentSummary);
    writetable(clusterContrasts, paths.clusterContrasts);
    writetable(clusterSummary, paths.clusterSummary);

    fprintf('Saved deltaBias permutation audit outputs:\n');
    fprintf('  %s\n', paths.mat);
    fprintf('  %s\n', paths.experimentContrasts);
    fprintf('  %s\n', paths.experimentSummary);
    fprintf('  %s\n', paths.clusterContrasts);
    fprintf('  %s\n', paths.clusterSummary);
end