function [figHandle, demoPlotData, outputFiles] = plotDemoOptostim( ...
    currentBlockStruct, imagingData, bitmapData, ...
    columnarProducts, demoProducts, blockID, currentSessID, ...
    mainPath, monkeyName, chamberWanted, saveFlag)
%PLOTDEMOOPTOSTIM  1x3 demonstration figure for a single optostim session.
%
% Panel A – PCA(90°-0°) background with ROI outline and Gaussian footprint
% Panel B – PCA(90°-0°) background with blue/red targeted-column overlays
% Panel C – Combined view (A + B)
%
% All panels share the same CLim and physical mm axes.
% Returns demoPlotData with all underlying numerical arrays for inspection.
%
% Usage:
%   [figHandle, demoPlotData, outputFiles] = plotDemoOptostim( ...
%       currentBlockStruct, imagingData, bitmapData, ...
%       columnarProducts, demoProducts, blockID, currentSessID, ...
%       mainPath, monkeyName, chamberWanted, saveFlag)

% ──────────────────────────────────────────────────────────────────
% 1. Orientation lookup (raises a descriptive error if 0° or 90° absent)
% ──────────────────────────────────────────────────────────────────
[idx0, idx90, bitmapIdx0, bitmapIdx90] = demoFindOrientationIndices( ...
    columnarProducts.orts, demoProducts.orts, currentSessID);

% ──────────────────────────────────────────────────────────────────
% 2. Extract PCA slices (reference-session space)
% ──────────────────────────────────────────────────────────────────
pca0  = columnarProducts.VERpca(:,:,idx0);
pca90 = columnarProducts.VERpca(:,:,idx90);

if all(isnan(pca0(:)))
    error('plotDemoOptostim:NaN0', 'Session %d: 0%s PCA map is all-NaN.', ...
        currentSessID, char(176));
end
if all(isnan(pca90(:)))
    error('plotDemoOptostim:NaN90', 'Session %d: 90%s PCA map is all-NaN.', ...
        currentSessID, char(176));
end

% ──────────────────────────────────────────────────────────────────
% 3. Warp PCA maps from reference space to current-session camera space
% ──────────────────────────────────────────────────────────────────
target0cam  = demoProducts.targetedColumnsCamspace(:,:,bitmapIdx0);
target90cam = demoProducts.targetedColumnsCamspace(:,:,bitmapIdx90);

[pca0Coreg, pca90Coreg, validMask] = demoWarpPca( ...
    pca0, pca90, bitmapData, blockID, target0cam);

pcaDifference = pca90Coreg - pca0Coreg;
pcaDifference(~validMask) = NaN;

% ──────────────────────────────────────────────────────────────────
% 4. Prepare binary masks
% ──────────────────────────────────────────────────────────────────
roiMask    = demoProducts.roiMaskCamspace;                                % logical
gaussMask0 = demoProducts.gaussianMaskCamspace(:,:,bitmapIdx0) > 0.5;    % logical
target0    = isfinite(target0cam)  & (target0cam  ~= 0);                  % logical
target90   = isfinite(target90cam) & (target90cam ~= 0);                  % logical

if ~any(roiMask(:))
    error('plotDemoOptostim:EmptyROI', 'Session %d: ROI mask is empty.', currentSessID);
end
if ~any(gaussMask0(:))
    warning('plotDemoOptostim:EmptyGauss', ...
        'Session %d: Gaussian footprint mask for 0%s is empty.', ...
        currentSessID, char(176));
end
if ~any(target0(:))
    error('plotDemoOptostim:EmptyTarget0', ...
        'Session %d: 0%s targeted-column bitmap is empty.', ...
        currentSessID, char(176));
end
if ~any(target90(:))
    error('plotDemoOptostim:EmptyTarget90', ...
        'Session %d: 90%s targeted-column bitmap is empty.', ...
        currentSessID, char(176));
end

% ──────────────────────────────────────────────────────────────────
% 5. Validate dimensions match
% ──────────────────────────────────────────────────────────────────
demoValidateDimensions(pcaDifference, roiMask, gaussMask0, target0, target90);

% ──────────────────────────────────────────────────────────────────
% 6. Physical mm coordinate vectors
% ──────────────────────────────────────────────────────────────────
[xMM, yMM, pixSize] = demoComputeMMCoords(imagingData, blockID, pcaDifference);

% ──────────────────────────────────────────────────────────────────
% 7. Symmetric color limits (calculated once, shared across all panels)
% ──────────────────────────────────────────────────────────────────
validPix = validMask & isfinite(pcaDifference);
if ~any(validPix(:))
    error('plotDemoOptostim:AllNaN', ...
        'Session %d: PCA difference map has no finite pixels.', currentSessID);
end
maxAbs = max(abs(pcaDifference(validPix)));
if ~isfinite(maxAbs) || maxAbs == 0
    error('plotDemoOptostim:ZeroRange', ...
        'Session %d: PCA difference map has zero dynamic range.', currentSessID);
end
cLim = [-maxAbs, maxAbs];

% ──────────────────────────────────────────────────────────────────
% 8. Check for column overlap
% ──────────────────────────────────────────────────────────────────
overlap  = target0 & target90;
nOverlap = sum(overlap(:));
if nOverlap > 0
    warning('plotDemoOptostim:ColumnOverlap', ...
        'Session %d: %d pixels overlap between 0%s and 90%s targeted columns. Shown in magenta.', ...
        currentSessID, nOverlap, char(176), char(176));
end

% ──────────────────────────────────────────────────────────────────
% 9. Diagnostics
% ──────────────────────────────────────────────────────────────────
fprintf('\n--- Demo optostim diagnostics ---\n');
fprintf('Session:                  %d\n',  currentSessID);
fprintf('PCA orientations found:   [%s]\n', num2str(columnarProducts.orts));
fprintf('Bitmap orientations found:[%s]\n', num2str(demoProducts.orts));
fprintf('Image dimensions:         %d x %d\n', size(pcaDifference,1), size(pcaDifference,2));
fprintf('Pixel size (mm/px):       %.6f\n', pixSize);
fprintf('ROI pixels:               %d\n',  sum(roiMask(:)));
fprintf('Gaussian footprint pixels:%d\n',  sum(gaussMask0(:)));
fprintf('0-deg target pixels:      %d\n',  sum(target0(:)));
fprintf('90-deg target pixels:     %d\n',  sum(target90(:)));
fprintf('Target overlap pixels:    %d\n',  nOverlap);
fprintf('PCA difference range:     [%.4g  %.4g]\n', ...
    min(pcaDifference(validPix)), max(pcaDifference(validPix)));
fprintf('---------------------------------\n\n');

% ──────────────────────────────────────────────────────────────────
% 10. Create figure
% ──────────────────────────────────────────────────────────────────
figHandle = figure( ...
    'Name',     sprintf('demo-optostim_sess%d', currentSessID), ...
    'Color',    'white', ...
    'Units',    'inches', ...
    'Position', [0.5  0.5  18  6]);

% Panel A
ax1 = subplot(1, 3, 1);
demoBuildPanelA(ax1, pcaDifference, xMM, yMM, cLim, pixSize, roiMask, gaussMask0);
text(ax1, -0.15, 1.04, 'A', 'Units', 'normalized', ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');

% Panel B
ax2 = subplot(1, 3, 2);
demoBuildPanelB(ax2, pcaDifference, xMM, yMM, cLim, target0, target90, nOverlap);
text(ax2, -0.15, 1.04, 'B', 'Units', 'normalized', ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');

% Panel C
ax3 = subplot(1, 3, 3);
demoBuildPanelC(ax3, pcaDifference, xMM, yMM, cLim, pixSize, ...
    roiMask, gaussMask0, target0, target90, nOverlap);
text(ax3, -0.15, 1.04, 'C', 'Units', 'normalized', ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');

% ──────────────────────────────────────────────────────────────────
% 11. Signed target map helper
% ──────────────────────────────────────────────────────────────────
signedTargetMap = nan(size(target0), 'double');
signedTargetMap(target0)  = -1;
signedTargetMap(target90) =  1;

% ──────────────────────────────────────────────────────────────────
% 12. Pack demoPlotData for caller inspection
% ──────────────────────────────────────────────────────────────────
demoPlotData.pca0Coreg      = pca0Coreg;
demoPlotData.pca90Coreg     = pca90Coreg;
demoPlotData.pcaDifference  = pcaDifference;
demoPlotData.roiMask        = roiMask;
demoPlotData.gaussianMask0  = gaussMask0;
demoPlotData.targetMask0    = target0;
demoPlotData.targetMask90   = target90;
demoPlotData.signedTargetMap = signedTargetMap;
demoPlotData.xMM            = xMM;
demoPlotData.yMM            = yMM;
demoPlotData.cLim           = cLim;

% ──────────────────────────────────────────────────────────────────
% 13. Save SVG and PNG
% ──────────────────────────────────────────────────────────────────
outputFiles = demoSaveOutputs(figHandle, mainPath, monkeyName, chamberWanted, ...
    currentSessID, currentBlockStruct, saveFlag);

end % plotDemoOptostim


% ==============================================================
%   HELPER SUBFUNCTIONS
% ==============================================================

function [idx0, idx90, bitmapIdx0, bitmapIdx90] = demoFindOrientationIndices( ...
        pcaOrts, bmpOrts, currentSessID)
% Find 0° and 90° in both orientation arrays; raise on missing/duplicate.

    idx0  = find(pcaOrts == 0,  1);
    idx90 = find(pcaOrts == 90, 1);

    if isempty(idx0)
        error('plotDemoOptostim:MissingPCA0', ...
            'Session %d: 0%s orientation absent from PCA data. Found: [%s]', ...
            currentSessID, char(176), num2str(pcaOrts));
    end
    if isempty(idx90)
        error('plotDemoOptostim:MissingPCA90', ...
            'Session %d: 90%s orientation absent from PCA data. Found: [%s]', ...
            currentSessID, char(176), num2str(pcaOrts));
    end
    if sum(pcaOrts == 0) > 1 || sum(pcaOrts == 90) > 1
        warning('plotDemoOptostim:DuplicatePCAOrts', ...
            'Session %d: duplicate orientations in PCA data; using first occurrence.', ...
            currentSessID);
    end

    bitmapIdx0  = find(bmpOrts == 0,  1);
    bitmapIdx90 = find(bmpOrts == 90, 1);

    if isempty(bitmapIdx0)
        error('plotDemoOptostim:MissingBitmap0', ...
            'Session %d: 0%s orientation absent from bitmap data. Found: [%s]', ...
            currentSessID, char(176), num2str(bmpOrts));
    end
    if isempty(bitmapIdx90)
        error('plotDemoOptostim:MissingBitmap90', ...
            'Session %d: 90%s orientation absent from bitmap data. Found: [%s]', ...
            currentSessID, char(176), num2str(bmpOrts));
    end
end


function [pca0Coreg, pca90Coreg, validMask] = demoWarpPca( ...
        pca0, pca90, bitmapData, blockID, refImage)
% Warp PCA maps from reference-session space into current-session camera space.
% Uses the same transform saved by coregisterBitmap2GreenImgV2.

    if ~isfield(bitmapData, 'transformParams') || ...
            blockID > numel(bitmapData.transformParams) || ...
            isempty(bitmapData.transformParams{blockID})
        error('plotDemoOptostim:MissingTransform', ...
            'bitmapData.transformParams{%d} is missing or empty.', blockID);
    end

    tf        = bitmapData.transformParams{blockID};
    outputRef = imref2d(size(refImage));

    pca0Coreg  = imwarp(pca0,  tf, 'OutputView', outputRef);
    pca90Coreg = imwarp(pca90, tf, 'OutputView', outputRef);

    % Pixels outside the warped field-of-view get interpolated zeros; mark them NaN.
    validRef   = ones(size(pca0), 'double');
    validCoreg = imwarp(validRef, tf, 'OutputView', outputRef);
    validMask  = validCoreg > 0.5;
end


function demoValidateDimensions(pcaDiff, roiMask, gaussMask, target0, target90)
% Verify all arrays share the same height and width.

    baseSize = size(pcaDiff);
    arrays   = {roiMask, gaussMask, target0, target90};
    names    = {'roiMask', 'gaussMask0', 'target0', 'target90'};

    for k = 1:numel(arrays)
        if size(arrays{k}, 1) ~= baseSize(1) || size(arrays{k}, 2) ~= baseSize(2)
            error('plotDemoOptostim:DimMismatch', ...
                'Dimension mismatch: pcaDifference is [%d %d] but %s is [%d %d].', ...
                baseSize(1), baseSize(2), names{k}, ...
                size(arrays{k}, 1), size(arrays{k}, 2));
        end
    end
end


function [xMM, yMM, pixSize] = demoComputeMMCoords(imagingData, blockID, refArray)
% Return physical mm coordinate vectors matching rows/columns of refArray.

    DEFAULT_MM_PER_PX = 8.22 / 512;   % addPix2MM convention

    pixSize = [];
    if isfield(imagingData, 'pixelsizemm') && ...
            numel(imagingData.pixelsizemm) >= blockID && ...
            isfinite(imagingData.pixelsizemm(blockID)) && ...
            imagingData.pixelsizemm(blockID) > 0
        pixSize = imagingData.pixelsizemm(blockID);
    end

    if isempty(pixSize)
        warning('plotDemoOptostim:NoPixelSize', ...
            'imagingData.pixelsizemm unavailable; using default %.6f mm/px (8.22 mm / 512 px).', ...
            DEFAULT_MM_PER_PX);
        pixSize = DEFAULT_MM_PER_PX;
    end

    [nRows, nCols] = size(refArray);
    xMM = (0:(nCols-1)) * pixSize;
    yMM = (0:(nRows-1)) * pixSize;
end


% ─── Panel builders ────────────────────────────────────────────────

function demoBuildPanelA(ax, pcaDiff, xMM, yMM, cLim, pixSize, roiMask, gaussMask0)
% Panel A: diverging PCA background + ROI outline + Gaussian footprint outline.

    demoPcaBackground(ax, pcaDiff, xMM, yMM, cLim);
    hold(ax, 'on');

    % ROI outline – black dotted
    demoOutline(ax, roiMask, xMM, yMM, pixSize, [0 0 0], ':',  1.5, 'ROI');
    % Gaussian footprint outline – orange dashed
    demoOutline(ax, gaussMask0, xMM, yMM, pixSize, [1 0.65 0], '--', 1.5, 'Gaussian footprint');

    demoFormatAxes(ax, xMM, yMM);

    cb = colorbar(ax);
    cb.Label.String = ['PCA response: 90' char(176) ' - 0' char(176)];
    cb.Label.FontSize  = 10;
    cb.Label.FontName  = 'Arial';

    title(ax, ['90' char(176) ' - 0' char(176) ' PCA response with ROI'], ...
        'FontWeight', 'normal', 'FontName', 'Arial', 'FontSize', 12);

    legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
        'FontSize', 9, 'FontName', 'Arial', 'Box', 'off');
end


function demoBuildPanelB(ax, pcaDiff, xMM, yMM, cLim, target0, target90, nOverlap)
% Panel B: diverging PCA background + blue/red column overlays.

    demoPcaBackground(ax, pcaDiff, xMM, yMM, cLim);
    hold(ax, 'on');

    demoRGBOverlay(ax, target0, target90, xMM, yMM, nOverlap);

    demoFormatAxes(ax, xMM, yMM);

    cb = colorbar(ax);
    cb.Label.String = ['PCA response: 90' char(176) ' - 0' char(176)];
    cb.Label.FontSize = 10;
    cb.Label.FontName = 'Arial';

    title(ax, ['Targeted 0' char(176) ' and 90' char(176) ' columns'], ...
        'FontWeight', 'normal', 'FontName', 'Arial', 'FontSize', 12);

    legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
        'FontSize', 9, 'FontName', 'Arial', 'Box', 'off');
end


function demoBuildPanelC(ax, pcaDiff, xMM, yMM, cLim, pixSize, ...
        roiMask, gaussMask0, target0, target90, nOverlap)
% Panel C: all overlays combined.

    demoPcaBackground(ax, pcaDiff, xMM, yMM, cLim);
    hold(ax, 'on');

    demoRGBOverlay(ax, target0, target90, xMM, yMM, nOverlap);
    demoOutline(ax, roiMask,   xMM, yMM, pixSize, [0 0 0],    ':',  1.5, 'ROI');
    demoOutline(ax, gaussMask0, xMM, yMM, pixSize, [1 0.65 0], '--', 1.5, 'Gaussian footprint');

    demoFormatAxes(ax, xMM, yMM);

    title(ax, 'Combined optostimulation targeting', ...
        'FontWeight', 'normal', 'FontName', 'Arial', 'FontSize', 12);

    legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
        'FontSize', 9, 'FontName', 'Arial', 'Box', 'off');
end


% ─── Primitive drawing helpers ──────────────────────────────────────

function demoPcaBackground(ax, pcaDiff, xMM, yMM, cLim)
% Render the diverging PCA-difference map as the background image.

    hImg = imagesc(ax, xMM, yMM, pcaDiff);
    hImg.HandleVisibility = 'off';   % keep out of auto-legend
    colormap(ax, redblue(256));
    set(ax, 'CLim', cLim);
end


function demoOutline(ax, mask, xMM, yMM, pixSize, color, lstyle, lwidth, label)
% Plot the binary-mask boundary as a line overlay using bwboundaries.
% Only the first boundary segment gets a legend entry.

    if ~any(mask(:))
        % Still add invisible proxy so legend entry exists (line not plotted)
        plot(ax, NaN, NaN, lstyle, 'Color', color, 'LineWidth', lwidth, ...
            'DisplayName', label);
        return;
    end

    bounds = bwboundaries(mask);  % {[row col]} pairs
    for k = 1:numel(bounds)
        b       = bounds{k};
        xCoords = (b(:,2) - 1) * pixSize;   % column index → mm
        yCoords = (b(:,1) - 1) * pixSize;   % row index    → mm
        if k == 1
            plot(ax, xCoords, yCoords, lstyle, ...
                'Color', color, 'LineWidth', lwidth, 'DisplayName', label);
        else
            plot(ax, xCoords, yCoords, lstyle, ...
                'Color', color, 'LineWidth', lwidth, 'HandleVisibility', 'off');
        end
    end
end


function demoRGBOverlay(ax, target0, target90, xMM, yMM, nOverlap)
% Overlay targeted columns as a semi-transparent RGB image.
%   0°  pixels → blue  [0 0 1]
%   90° pixels → red   [1 0 0]
%   overlap    → magenta [1 0 1]

    ALPHA    = 0.72;
    COL0     = [0 0 1];
    COL90    = [1 0 0];
    COLOVLP  = [1 0 1];

    [nR, nC] = size(target0);
    rgbIm    = zeros(nR, nC, 3, 'double');
    alphaIm  = zeros(nR, nC, 'double');

    % 0° – blue
    for ch = 1:3
        plane = rgbIm(:,:,ch);
        plane(target0) = COL0(ch);
        rgbIm(:,:,ch) = plane;
    end
    alphaIm(target0) = ALPHA;

    % 90° – red (non-overlap pixels only; overlap handled below)
    noOvlp90 = target90 & ~(target0 & target90);
    for ch = 1:3
        plane = rgbIm(:,:,ch);
        plane(noOvlp90) = COL90(ch);
        rgbIm(:,:,ch) = plane;
    end
    alphaIm(target90) = ALPHA;

    % Overlap – magenta
    if nOverlap > 0
        ovlp = target0 & target90;
        for ch = 1:3
            plane = rgbIm(:,:,ch);
            plane(ovlp) = COLOVLP(ch);
            rgbIm(:,:,ch) = plane;
        end
        alphaIm(ovlp) = 0.90;
    end

    hOvlp = image(ax, xMM, yMM, rgbIm);
    hOvlp.AlphaData        = alphaIm;
    hOvlp.HandleVisibility = 'off';

    % Proxy patches for the legend (invisible in the plot)
    patch(ax, NaN, NaN, COL0, 'FaceAlpha', ALPHA, 'EdgeColor', 'none', ...
        'DisplayName', ['0' char(176) ' targeted columns']);
    patch(ax, NaN, NaN, COL90, 'FaceAlpha', ALPHA, 'EdgeColor', 'none', ...
        'DisplayName', ['90' char(176) ' targeted columns']);
    if nOverlap > 0
        patch(ax, NaN, NaN, COLOVLP, 'FaceAlpha', 0.90, 'EdgeColor', 'none', ...
            'DisplayName', sprintf('Overlap (%d px)', nOverlap));
    end
end


function demoFormatAxes(ax, xMM, yMM)
% Apply consistent axes formatting to every demo panel.

    axis(ax, 'image');
    set(ax, 'YDir',    'reverse');
    set(ax, 'TickDir', 'out');
    set(ax, 'FontName', 'Arial', 'FontSize', 11);
    set(ax, 'Box', 'off');
    set(ax, 'XGrid', 'off', 'YGrid', 'off');
    xlim(ax, [xMM(1)  xMM(end)]);
    ylim(ax, [yMM(1)  yMM(end)]);
    xlabel(ax, 'X (mm)', 'FontName', 'Arial');
    ylabel(ax, 'Y (mm)', 'FontName', 'Arial');
end


% ─── Save helper ────────────────────────────────────────────────────

function outputFiles = demoSaveOutputs(figHandle, mainPath, monkeyName, ...
        chamberWanted, currentSessID, currentBlockStruct, saveFlag)
% Save SVG and PNG (300 dpi) when saveFlag == 1.

    outputFiles = struct('svg', '', 'png', '');

    if saveFlag ~= 1
        return;
    end

    % Output directory
    outDir = fullfile(mainPath, monkeyName, 'Meta', 'demo-optostim');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    % Filesystem-safe filename components
    dateStr = regexprep(char(currentBlockStruct.date), '[^A-Za-z0-9_-]', '');
    runStr  = regexprep(num2str(currentBlockStruct.run), '[^A-Za-z0-9_-]', '');
    chamStr = regexprep(chamberWanted, '[^A-Za-z0-9_-]', '');

    baseName = sprintf('demo-optostim_%s_%s_sess%d_%sR%s', ...
        monkeyName, chamStr, currentSessID, dateStr, runStr);

    svgFile = fullfile(outDir, [baseName '.svg']);
    pngFile = fullfile(outDir, [baseName '.png']);

    set(figHandle, 'PaperPositionMode', 'auto');

    print(figHandle, svgFile, '-dsvg');
    print(figHandle, pngFile, '-dpng', '-r300');

    fprintf('Saved SVG: %s\n', svgFile);
    fprintf('Saved PNG: %s\n', pngFile);

    outputFiles.svg = svgFile;
    outputFiles.png = pngFile;
end
