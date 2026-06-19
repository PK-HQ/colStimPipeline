function reportState = stageReportPDFPage(reportState, figHandle)
% Export one live figure to an ordered temporary PDF page.

    if isempty(figHandle) || ~isgraphics(figHandle)
        error(['Plot producer did not return a live graphics handle. ' ...
            'class=%s, size=%s'], class(figHandle), mat2str(size(figHandle)));
    end
    if ~isfield(reportState, 'tempDir') || ...
            exist(reportState.tempDir, 'dir') ~= 7
        error('Report temporary directory is missing.');
    end

    reportState.pageIdx = reportState.pageIdx + 1;
    pageFilename = fullfile(reportState.tempDir, ...
        sprintf('page_%04d.pdf', reportState.pageIdx));
    [~, ~, pageExt] = fileparts(pageFilename);
    assert(strcmpi(pageExt, '.pdf'), ...
        'Temporary page filename must end in .pdf: %s', pageFilename);

    drawnow;
    exportgraphics(figHandle, pageFilename, ...
        'ContentType', 'image', ...
        'Resolution', reportState.resolution, ...
        'BackgroundColor', 'white');
    assertPdfExists(pageFilename, 'temporary page');

    reportState.pageFiles{reportState.pageIdx} = pageFilename;
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
