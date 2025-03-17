function scaledParams = scaleParameters(params, lb, ub)
% SCALEPARAMETERS - Scale parameters to [0,1] range for optimization
%
% Inputs:
%   params - Original parameters in their natural units
%   lb, ub - Lower and upper bounds for parameters
%
% Outputs:
%   scaledParams - Parameters scaled to [0,1] range

    % Handle different parameter scale ranges
    scaledParams = zeros(size(params));
    
    for i = 1:length(params)
        if i >= 5 && i <= 7  % g0, n, rmx - log scale
            % Convert to log scale for wide-range parameters
            logLB = log(lb(i));
            logUB = log(ub(i));
            logParam = log(params(i));
            
            % Scale to [0,1]
            scaledParams(i) = (logParam - logLB) / (logUB - logLB);
        else  % linear scale for other parameters
            % Linear scaling to [0,1]
            scaledParams(i) = (params(i) - lb(i)) / (ub(i) - lb(i));
        end
        
        % Ensure in [0,1] range (handle numerical precision issues)
        scaledParams(i) = max(0, min(1, scaledParams(i)));
    end
end