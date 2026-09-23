function reportState = stageReportPDFPage(reportState, figHandle)
% Append one live figure to the single local report PDF.

    pageIndex = reportState.pageIdx + 1;
    if ~isscalar(figHandle) || ~isgraphics(figHandle, 'figure')
        error('stageReportPDFPage:InvalidFigure', ...
            'PDF page %d requires one live figure handle.', pageIndex);
    end
    if exist(reportState.tempDir, 'dir') ~= 7
        error('stageReportPDFPage:MissingDirectory', ...
            'Report temporary directory is missing: %s', reportState.tempDir);
    end
    if ~isempty(reportState.expectedPages) && pageIndex > reportState.expectedPages
        error('stageReportPDFPage:TooManyPages', 'More figures than requested.');
    end
    beforeBytes = 0;
    if reportState.pageIdx > 0
        info = dir(reportState.localPdf);
        if isempty(info) || info.bytes ~= reportState.pageBytes(end, 2)
            error('stageReportPDFPage:StaleState', ...
                'Local PDF changed or disappeared: %s', reportState.localPdf);
        end
        beforeBytes = info.bytes;
    elseif exist(reportState.localPdf, 'file') == 2
        error('stageReportPDFPage:StaleState', ...
            'Untracked PDF already exists: %s', reportState.localPdf);
    end

    % An interrupted export invalidates the assembly, even if the file survived.
    pendingFile = fullfile(reportState.tempDir, 'export_pending');
    if exist(pendingFile, 'file') == 2
        error('stageReportPDFPage:InterruptedExport', ...
            'An earlier export did not complete: %s', reportState.tempDir);
    end
    fid = fopen(pendingFile, 'w');
    if fid < 0
        error('stageReportPDFPage:PendingMarker', ...
            'Cannot write export marker: %s', pendingFile);
    end
    fclose(fid);
    try
        drawnow;
        if ~isgraphics(figHandle, 'figure')
            error('stageReportPDFPage:ClosedFigure', 'Figure closed during drawnow.');
        end
        exportgraphics(figHandle, reportState.localPdf, ...
            'ContentType', 'image', 'Resolution', reportState.resolution, ...
            'BackgroundColor', 'white', 'Append', pageIndex > 1);
        info = dir(reportState.localPdf);
        if isempty(info) || info.bytes <= beforeBytes
            error('stageReportPDFPage:NoGrowth', ...
                'PDF missing, empty, or did not grow after page %d.', pageIndex);
        end
    catch exportError
        error('stageReportPDFPage:ExportFailed', ...
            'PDF page %d failed. Local files preserved in %s. Cause: %s', ...
            pageIndex, reportState.tempDir, exportError.message);
    end
    delete(pendingFile);
    reportState.pageIdx = pageIndex;
    reportState.pageBytes(pageIndex, :) = [beforeBytes info.bytes];
    if isempty(reportState.expectedPages)
        totalText = '?';
    else
        totalText = num2str(reportState.expectedPages);
    end
    fprintf('PDF page %d/%s saved | bytes %d -> %d\n', ...
        pageIndex, totalText, beforeBytes, info.bytes);
end
