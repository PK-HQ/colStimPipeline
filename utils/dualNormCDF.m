% --- Curve fitting for data ---
function yFit = dualNormCDF(param,x)
% Version: Feb 2023
% Psychometric Function for contrast, 2 asymmetrical curves for x<mu and x>mu

muLower = param(1);
sigmaLower = param(2);
muUpper = param(3);
sigmaUpper = param(4);
beta = param(5);
lapseLower = param(6);
lapseUpper = param(7);

% Fit 2 parts (<0 and >0), combine
xLower=x(x<0);%mu);
xUpper=x(x>0);%mu);
yFitLower=((beta) * normcdf(xLower,-muLower,sigmaLower))+lapseLower;
yFitUpper=(((1-beta-lapseLower-lapseUpper) * normcdf(xUpper,muUpper,sigmaUpper))) + (beta+lapseLower);
% Combine to get single curve spanning x-range
yFit=[yFitLower, yFitUpper];

% clipping, carryover from Giac's script
yFit(yFit<0.0001) = 0.0001;
yFit(yFit> 0.9999) = 0.9999;

end