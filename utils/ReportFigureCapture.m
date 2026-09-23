classdef ReportFigureCapture < handle
% Record figure creation events without changing root defaults or figure order.
    properties (Access = private)
        Figures = gobjects(0, 1)
        Listener
        Finished = false
    end
    methods
        function obj = ReportFigureCapture()
            obj.Listener = addlistener(groot, 'ObjectChildAdded', ...
                @(~, event) obj.recordFigure(event.Child));
        end
        function figures = finish(obj)
            if obj.Finished
                error('ReportFigureCapture:AlreadyFinished', ...
                    'This figure capture has already ended.');
            end
            delete(obj.Listener);
            obj.Finished = true;
            figures = obj.Figures(isgraphics(obj.Figures, 'figure'));
        end
        function delete(obj)
            if ~isempty(obj.Listener)
                delete(obj.Listener);
            end
        end
    end
    methods (Access = private)
        function recordFigure(obj, child)
            if isgraphics(child, 'figure')
                obj.Figures(end + 1, 1) = child;
            end
        end
    end
end
