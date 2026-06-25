function [lb, ub] = getWeibullSignedBX0PanelBounds()
%GETWEIBULLSIGNEDBX0PANELBOUNDS Bounds for the 10-parameter conditional panel model.
%
% beta > 1 gives a continuous first derivative at onset;
% beta >= 2 prevents divergent onset curvature in the flat-then-rising
% conditional panel representation.
%
% The primary joint model bounds (11 parameters) are in getWeibullSignedBX0InitParams.
    [~, lbFull, ubFull] = getWeibullSignedBX0InitParams();
    lb = lbFull(1:10);
    ub = ubFull(1:10);
    lb([3, 6, 9]) = 2.0;
end
