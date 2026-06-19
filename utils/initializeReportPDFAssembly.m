function reportState = initializeReportPDFAssembly(filename, monkeyName, opts)
% Initialize staged PDF assembly in a local temporary directory.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'resolution') || isempty(opts.resolution)
        opts.resolution = 150;
    end

    reportState = struct();
    reportState.outputFilename = resolveReportFilename(filename, monkeyName);
    reportState.tempDir = fullfile(tempdir, ...
        ['colStimReport_' char(java.util.UUID.randomUUID)]);
    mkdir(reportState.tempDir);
    reportState.pageFiles = {};
    reportState.pageIdx = 0;
    reportState.resolution = opts.resolution;
end

function outputFilename = resolveReportFilename(filename, monkeyName)
    if endsWith(filename, '.pdf', 'IgnoreCase', true)
        outputFilename = filename;
    else
        if ispc
            mainPath = 'Y:/';
        elseif contains(getenv('HOSTNAME'), 'psy.utexas.edu')
            mainPath = '/eslab/data/';
        else
            error('Unable to determine the output root for PDF saving.');
        end
        outputFilename = [mainPath monkeyName '/Meta/' filename '.pdf'];
    end

    [~, ~, outputExt] = fileparts(outputFilename);
    assert(strcmpi(outputExt, '.pdf'), ...
        'Final report filename must end in .pdf: %s', outputFilename);
end
