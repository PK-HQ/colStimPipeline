function [initialParams, lb, ub] = getWeibullInitParams()
     initialParams = [ ...
        0.10  0.50  15  3 , ...   % baseline      (A_bl  B_bl  α_bl  β_bl)
        0.00  0.00   0  0 , ...   % Δ-congruent   (ΔA    ΔB    Δα   Δβ)
        0.00  0.00   0  0 ];      % Δ-incongruent (ΔA    ΔB*   Δα   Δβ)
        % *ΔB_inc is a dummy slot; it is fixed to 0 by the bounds below.
    lb = [ ...
        0.00  0.50  10  1 , ...    % baseline  (B_bl is fixed at 0.5)
       -0.20 -0.30 -10 -4 , ...    % Δ-congruent  (allow up/down shifts)
       -0.20  0.00 -10 -4 ];       % Δ-incongruent (ΔB_inc fixed at 0)
    
    ub = [ ...
        0.20  0.50  50  8 , ...    % baseline  (B_bl fixed at 0.5)
        0.20  0.30  10  4 , ...    % Δ-congruent
        0.20  0.00  10  4 ];       % Δ-incongruent (ΔB_inc fixed at 0)
end

        