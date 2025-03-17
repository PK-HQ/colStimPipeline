function cax2(varargin)
% varargin{1} = modifier (optional)
% varargin{2} = img (optional)

    if nargin == 0
        modifier = 0.8;
        caxisLimits = caxis; % Get current caxis limits
    elseif nargin == 1
        modifier = varargin{1};
        caxisLimits = caxis;
    elseif nargin >= 2
        modifier = varargin{1};
        img = varargin{2};
        caxisLimits = [min(img(:), [], 'omitnan'), max(img(:), [], 'omitnan')];
    end

    if isempty(caxisLimits) || (caxisLimits(2)<=caxisLimits(1))  % Handle cases where caxis is empty or invalid
        return; % Or you could set some default caxis limits here
    end

    if round(caxisLimits(1), 2, 'significant') < 0 && caxisLimits(2) > 0  % +/- range
        caxisLimitsNew = max(abs(caxisLimits));
        caxis([-caxisLimitsNew, caxisLimitsNew] * modifier);
    elseif round(caxisLimits(1), 2, 'significant') >= 0 && caxisLimits(2) > 0  % + range
        caxis([0, caxisLimits(2)] * modifier);
    elseif caxisLimits(1) < 0 && round(caxisLimits(2), 2, 'significant') <= 0  % - range
        caxis([caxisLimits(1), 0] * modifier);
    % No else needed. If none of the above, we keep the existing caxis limits.
    end

    % Add a colorbar (if one doesn't already exist)
    if isempty(colorbar('peer', gca))
      colorbar;
    end
    
    if nargin >=2 %resize colorbar if image provided, use limits determined by image min/max
      colorbarHandle=colorbar;
      colorbarHandle.Limits = [min(img(:), [], 'omitnan'), max(img(:), [], 'omitnan')];
    end
end

    