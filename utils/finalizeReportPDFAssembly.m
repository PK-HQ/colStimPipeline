function outputFilename = finalizeReportPDFAssembly(reportState)
% Verify the completed local PDF, then copy it once to the destination.

    if reportState.pageIdx < 1 || ...
            size(reportState.pageBytes, 1) ~= reportState.pageIdx || ...
            any(reportState.pageBytes(:, 2) <= reportState.pageBytes(:, 1)) || ...
            exist(fullfile(reportState.tempDir, 'export_pending'), 'file') == 2
        error('finalizeReportPDFAssembly:IncompleteReport', ...
            'No complete report is available. Local files: %s', reportState.tempDir);
    end
    requestedPages = reportState.expectedPages;
    if isempty(requestedPages)
        requestedPages = reportState.pageIdx;
    end
    if reportState.pageIdx ~= requestedPages
        error('finalizeReportPDFAssembly:PageCountMismatch', ...
            'Requested %d figures but saved %d pages. Local PDF: %s', ...
            requestedPages, reportState.pageIdx, reportState.localPdf);
    end
    if ~isfield(reportState, 'outputFilename') || ...
            isempty(reportState.outputFilename)
        error('Report destination filename is missing.');
    end

    localFinalPdf = reportState.localPdf;
    assertPdfExists(localFinalPdf, 'local final report');
    localInfo = dir(localFinalPdf);
    if localInfo.bytes ~= reportState.pageBytes(end, 2)
        error('finalizeReportPDFAssembly:ChangedPDF', ...
            'Local PDF changed after its last successful export: %s', localFinalPdf);
    end
    fprintf(['PDF finalization | figures requested=%d | pages written=%d\n' ...
        'Local PDF: %s\nDestination: %s\nFinal bytes: %d\n'], ...
        requestedPages, reportState.pageIdx, localFinalPdf, ...
        reportState.outputFilename, localInfo.bytes);

    outputFilename = reportState.outputFilename;
    outputFolder = fileparts(outputFilename);
    try
        if ~isempty(outputFolder) && exist(outputFolder, 'dir') ~= 7
            mkdir(outputFolder);
        end
        [copyOK, copyMessage] = copyfile(localFinalPdf, outputFilename, 'f');
        if ~copyOK
            error('finalizeReportPDFAssembly:CopyFailed', '%s', copyMessage);
        end
        assertPdfExists(outputFilename, 'destination report');
        destinationInfo = dir(outputFilename);
        if destinationInfo.bytes ~= localInfo.bytes
            error('finalizeReportPDFAssembly:CopySizeMismatch', ...
                'Destination byte size differs from the completed local PDF.');
        end
    catch copyError
        fprintf(2, ['Completed local report preserved at: %s\n' ...
            'Temporary report directory preserved at: %s\n'], ...
            localFinalPdf, reportState.tempDir);
        error('finalizeReportPDFAssembly:DestinationLocked', ...
            ['Unable to replace destination PDF. Close the destination ' ...
            'PDF if it is open and retry. copyfile message: %s'], ...
            copyError.message);
    end

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
