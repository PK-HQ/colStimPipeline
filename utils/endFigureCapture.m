function figHandles = endFigureCapture(capture)
% Return still-live figures in creation order, then stop listening.
    if ~isa(capture, 'ReportFigureCapture') || ~isscalar(capture) || ~isvalid(capture)
        error('endFigureCapture:InvalidCapture', 'Supply a live figure capture.');
    end
    figHandles = capture.finish();
end
