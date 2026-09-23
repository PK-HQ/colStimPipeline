function capture = beginFigureCapture()
% Start recording newly created figures. Existing figures are not captured.
% Keep capture alive until endFigureCapture; deleting it cancels the listener.
    capture = ReportFigureCapture();
end
