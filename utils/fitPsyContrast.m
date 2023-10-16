function [mu,sigma,beta,logLikelihood] = fitPsyContrast(x,y)
%% README
% Fits psychometric function for contrast, 2 asymmetrical curves for x<mu and x>mu. Maximum Likelihood fit. Created
% by PK (2023), modified from Giac's script.

% --- Input ---
%       x = contrasts, signed according to orientation (0=-ve, 90=+ve)
%       y = % vertical reports

% --- Output ---
%       beta = vertical bias, scaling factor that modulates the height of normal cdf
%       mu = midpoint of cdf
%       sigma = spread of cdf
%       logLikelihood = goodness of fit

%% Fit data, find optimal parameter values
% --- Initial values, lower and upper bounds --- 
mu = 10; %is signed during fitting
sigma = 5;
beta = 0.5;
initialParams=[mu          sigma      beta];
paramLB      =[5           1             0];
paramUB      =[20          10            1];

% --- Run parameter optimization --- 
opts = optimset('fminsearch');
opts.MaxFunEvals = 100;
opts.MaxIter = 100;
[param]=fminsearchbnd(@(param) costFunction(param,x,y),...
    initialParams,...
    paramLB,...
    paramUB);

%  --- Extract param values and likelihood of best fit --- 
mu = param(1);
sigma = param(2);
beta = param(3);
logLikelihood  = costFunction([mu sigma beta],x,y);
end


%% Cost function: to calculate likelihood of fit vs actual data
function cost = costFunction(param,x,y)
%sort y and inverse of y
yInverse=100-y;

% Fit curve to data
yFit = dualNormCDF([param(1),param(2),param(3)],x);

% calculate probability for each data point
ct = 0;
for i =1:numel(yFit)
    if yFit(i) == 0 | yFit(i) == 1
        err
    else
        ct = ct+1;
        LLcont(ct) =  y(:,i).*log(yFit(i)) + yInverse(:,i).*log(1-yFit(i))  ;
    end
end
LL =LLcont(:) ; % log likelihood
cost = -2*sum(LL) ;
end

%% Curve fitting
function yFit = dualNormCDF(param,x)

% Psychometric Function for contrast, 2 asymmetrical curves for x<mu and x>mu
mu = param(1);
sigma = param(2);
beta = param(3);

% Fit 2 parts (<0 and >0), combine
xLower=x(x<0);%mu);
xUpper=x(x>0);%mu);
yFitLower=((beta) * normcdf(xLower,-mu,sigma));
yFitUpper=(((1-beta) * normcdf(xUpper,mu,sigma))) + beta;

% Combine to get single curve spanning x-range
yFit=[yFitLower, yFitUpper];

% Clipping, carryover from Giac's script
yFit(yFit<0.0001) = 0.0001;
yFit(yFit> 0.9999) = 0.9999;

end