function outputFilename = saveCompressedPDFPage(filename, monkeyName, figHandle, appendPage, opts)
% Append a figure to a compact multipage PDF without temporary page files.

    if nargin < 5 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'resolution') || isempty(opts.resolution)
        opts.resolution = 150;
    end
    validateattributes(opts.resolution, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});

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

    outputFolder = fileparts(outputFilename);
    if exist(outputFolder, 'dir') ~= 7
        mkdir(outputFolder);
    end
    if ~appendPage && exist(outputFilename, 'file') == 2
        delete(outputFilename);
    end

    drawnow;
    exportgraphics(figHandle, outputFilename, ...
        'ContentType', 'image', ...
        'Resolution', opts.resolution, ...
        'BackgroundColor', 'white', ...
        'Append', logical(appendPage));
end
