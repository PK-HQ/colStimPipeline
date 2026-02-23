% ---- helper: returns [ledPercent , ndFilter]
function [ledPrct, ndFilt] = parseLEDstring(str)
    % capture the numbers after L and after OD
    tk = regexp(str,'L(\d+)OD(\d+)','tokens','once');
    if isempty(tk)
        error('String "%s" is not in the expected LxxODyy format.',str);
    end
    ledPrct = str2double(tk{1});      % 30   from L30OD16
    ndFilt  = str2double(tk{2})/10;   % 1.6  from OD16  (remove /10 if you prefer "16")
end
