function outputFilename = finalizeReportPDFAssembly(reportState)
% Merge staged PDF pages once, then copy the complete report to destination.

    if ~isfield(reportState, 'pageFiles') || isempty(reportState.pageFiles)
        error('No staged PDF pages are available to merge.');
    end
    if ~isfield(reportState, 'outputFilename') || ...
            isempty(reportState.outputFilename)
        error('Report destination filename is missing.');
    end

    localFinalPdf = fullfile(reportState.tempDir, 'complete_report.pdf');
    append_pdfs(localFinalPdf, reportState.pageFiles{:});
    assertPdfExists(localFinalPdf, 'local final report');

    outputFilename = reportState.outputFilename;
    outputFolder = fileparts(outputFilename);
    if exist(outputFolder, 'dir') ~= 7
        mkdir(outputFolder);
    end

    [copyOK, copyMessage] = copyfile(localFinalPdf, outputFilename, 'f');
    if ~copyOK
        fprintf(2, ['Completed local report preserved at: %s\n' ...
            'Temporary report directory preserved at: %s\n'], ...
            localFinalPdf, reportState.tempDir);
        error('finalizeReportPDFAssembly:DestinationLocked', ...
            ['Unable to replace destination PDF. Close the destination ' ...
            'PDF if it is open and retry. copyfile message: %s'], ...
            copyMessage);
    end
    assertPdfExists(outputFilename, 'network final report');

    try
        rmdir(reportState.tempDir, 's');
    catch cleanupError
        warning('finalizeReportPDFAssembly:CleanupFailed', ...
            'Could not delete temporary report directory %s: %s', ...
            reportState.tempDir, cleanupError.message);
    end
end

function assertPdfExists(pdfFile, description)
    if exist(pdfFile, 'file') ~= 2
        error('Expected %s PDF does not exist: %s', description, pdfFile);
    end

    fileInfo = dir(pdfFile);
    if isempty(fileInfo) || fileInfo.bytes <= 0
        error('Expected %s PDF is empty: %s', description, pdfFile);
    end
end
