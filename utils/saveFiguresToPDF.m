function [outputPath, reportState] = saveFiguresToPDF(figHandles, outputFilename, opts)
% Save an ordered vector or cell array of live figures, one figure per page.
% Does not close figures. opts.resolution defaults to 150 dpi.
% The optional second output contains the successful count and byte audit.

    if nargin < 3
        opts = struct();
    end
    if isempty(figHandles) || ~isvector(figHandles)
        error('saveFiguresToPDF:InvalidCollection', ...
            'Supply a nonempty ordered vector or cell array of figures.');
    end
    for idx = 1:numel(figHandles)
        if iscell(figHandles)
            fig = figHandles{idx};
        else
            fig = figHandles(idx);
        end
        if ~isscalar(fig) || ~isgraphics(fig, 'figure')
            error('saveFiguresToPDF:InvalidFigure', ...
                'Requested page %d is not a live figure.', idx);
        end
    end
    opts.expectedPages = numel(figHandles);
    reportState = initializeReportPDFAssembly(outputFilename, '', opts);
    for idx = 1:numel(figHandles)
        if iscell(figHandles)
            fig = figHandles{idx};
        else
            fig = figHandles(idx);
        end
        reportState = stageReportPDFPage(reportState, fig);
    end
    outputPath = finalizeReportPDFAssembly(reportState);
end
