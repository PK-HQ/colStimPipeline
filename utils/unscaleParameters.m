function origParams = unscaleParameters(scaledParams, lb, ub)
% UNSCALEPARAMETERS - Convert scaled [0,1] parameters back to original scale
%
% Inputs:
%   scaledParams - Parameters in [0,1] range
%   lb, ub - Lower and upper bounds for parameters
%
% Outputs:
%   origParams - Parameters in their original units

    origParams = zeros(size(scaledParams));
    
    for i = 1:length(scaledParams)
        if i >= 5 && i <= 7  % g0, n, rmx - log scale
            % Convert from log scale for wide-range parameters
            logLB = log(lb(i));
            logUB = log(ub(i));
            logParam = logLB + scaledParams(i) * (logUB - logLB);
            
            % Back to original scale
            origParams(i) = exp(logParam);
        else  % linear scale for other parameters
            % Linear scaling from [0,1]
            origParams(i) = lb(i) + scaledParams(i) * (ub(i) - lb(i));
        end
    end
end