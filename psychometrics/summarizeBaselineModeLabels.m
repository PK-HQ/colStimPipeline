function textValue = summarizeBaselineModeLabels(modes)
% Summarize combined/separate baseline labels for plot metadata text.

    modes = string(modes(:));
    modes = modes(strlength(modes) > 0);
    if isempty(modes)
        textValue = 'n/a';
        return;
    end

    uniqueModes = unique(modes, 'stable');
    if numel(uniqueModes) == 1
        textValue = char(uniqueModes);
        return;
    end

    parts = strings(numel(uniqueModes), 1);
    for idx = 1:numel(uniqueModes)
        parts(idx) = sprintf('%s n=%d', ...
            uniqueModes(idx), sum(modes == uniqueModes(idx)));
    end
    textValue = char(strjoin(parts, ', '));
end
