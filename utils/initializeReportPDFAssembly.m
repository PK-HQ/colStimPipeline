function reportState = initializeReportPDFAssembly(filename, monkeyName, opts)
% Initialize staged PDF assembly in a local temporary directory.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'resolution') || isempty(opts.resolution)
        opts.resolution = 150;
    end

    validateattributes(opts.resolution, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    if nargin < 2
        monkeyName = '';
    end
    reportState = struct();
    reportState.outputFilename = resolveReportFilename(filename, monkeyName);
    reportState.tempDir = fullfile(tempdir, ...
        ['colStimReport_' char(java.util.UUID.randomUUID)]);
    mkdir(reportState.tempDir);
    reportState.localPdf = fullfile(reportState.tempDir, 'complete_report.pdf');
    reportState.pageBytes = zeros(0, 2);
    reportState.pageIdx = 0;
    reportState.expectedPages = [];
    if isfield(opts, 'expectedPages')
        validateattributes(opts.expectedPages, {'numeric'}, ...
            {'scalar', 'integer', 'positive', 'finite'});
        reportState.expectedPages = opts.expectedPages;
    end
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
