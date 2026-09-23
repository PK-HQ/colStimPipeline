function outputDir = smokeTestSaveFiguresToPDF(pdfinfoExecutable)
% Focused synthetic test; never loads analysis data or changes run settings.
% Optional pdfinfo executable verifies actual PDF page counts.
    if nargin < 1
        pdfinfoExecutable = '';
    end
    outputDir = tempname;
    mkdir(outputDir);
    fprintf('PDF smoke-test output: %s\n', outputDir);
    for n = [1 3 7]
        figures = gobjects(n, 1);
        for idx = 1:n
            figures(idx) = figure('Visible', 'off', 'Color', 'white', ...
                'Position', [100 100 480 320]);
            ax = axes('Parent', figures(idx));
            text(ax, 0.5, 0.5, sprintf('Page %d of %d', idx, n), ...
                'HorizontalAlignment', 'center', 'FontSize', 28);
            axis(ax, [0 1 0 1]);
            axis(ax, 'off');
        end
        cleanup = onCleanup(@() delete(figures(isgraphics(figures, 'figure'))));
        destination = fullfile(outputDir, sprintf('pages_%d.pdf', n));
        [path, state] = saveFiguresToPDF(figures, destination);
        assert(state.pageIdx == n && state.expectedPages == n);
        assert(size(state.pageBytes, 1) == n);
        assert(all(state.pageBytes(:, 2) > state.pageBytes(:, 1)));
        assert(state.pageBytes(1, 1) == 0);
        assert(isequal(state.pageBytes(2:end, 1), state.pageBytes(1:end-1, 2)));
        info = dir(path);
        assert(~isempty(info) && info.bytes == state.pageBytes(end, 2));
        assert(all(isgraphics(figures, 'figure')));
        if ~isempty(pdfinfoExecutable)
            [status, details] = system(sprintf('"%s" "%s"', pdfinfoExecutable, path));
            assert(status == 0, 'pdfinfo failed: %s', details);
            count = regexp(details, '(?m)^Pages:\s+(\d+)', 'tokens', 'once');
            assert(~isempty(count) && str2double(count{1}) == n, ...
                'Actual PDF page count differs from requested count.');
        end
        clear cleanup;
    end

    existing = figure('Visible', 'off');
    existingCleanup = onCleanup(@() delete(existing(isgraphics(existing))));
    capture = beginFigureCapture();
    captureCleanup = onCleanup(@() delete(capture));
    first = figure('Visible', 'off', 'HandleVisibility', 'off');
    firstCleanup = onCleanup(@() delete(first(isgraphics(first))));
    discarded = figure('Visible', 'off');
    delete(discarded);
    last = figure('Visible', 'off');
    lastCleanup = onCleanup(@() delete(last(isgraphics(last))));
    set(groot, 'CurrentFigure', first);
    captured = endFigureCapture(capture);
    assert(isequal(captured, [first; last]), ...
        'Capture must exclude old/deleted figures and preserve creation order.');

    invalidPath = fullfile(outputDir, 'invalid.pdf');
    try
        saveFiguresToPDF({first, discarded}, invalidPath);
        error('smokeTestSaveFiguresToPDF:MissingError', 'Invalid figure was accepted.');
    catch err
        assert(strcmp(err.identifier, 'saveFiguresToPDF:InvalidFigure'), '%s', err.message);
    end
    assert(exist(invalidPath, 'file') ~= 2);
    fprintf('PASS: 1/3/7 pages, byte growth, final page, capture, invalid handles.\n');
end
