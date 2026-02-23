%% 1. Load data
% Some 2-D FFTed or integrated array of functional imaging data
load('Y:\Pepper\Pepper20250708\run1\M32D20250708R1StabDFFTAmpS004E023PF0200.mat')
gaussImg = DataCond(:,:,2) - DataCond(:,:,1);               % Gaussian – blank

% The "maximal" response over space; usually taken from flashed large
% grating-evoked responses
load('Y:\Pepper\Pepper20250708\run0\M32D20250708R0StabDFFTAmpS004E044PF0200.mat')
refMap = DataCond(:,:,3:end) - mean(DataCond(:,:,1:2),3);  % 12‑ori grating – blank

% Scale the map to max out at 1
% You don't actually have to do this; if you don't then your fitted 
% Gaussian's amplitude will just be 1 or something
refMap = refMap ./ max(refMap(:));

% Define coordinate space
coord.X = 1:size(gaussImg, 1);
coord.Y = 1:size(gaussImg, 2);
coord.expression = refMap;

% Define initial 2D Gaussian parameters
PI  = [0.5,                     ... % amplitude
       ceil(size(gaussImg, 1) / 2), ... % center X coordinate
       ceil(size(gaussImg, 2) / 2), ... % center Y coordinate
       0,                       ... % DC component
       0,                       ... % orientation
       size(gaussImg, 1) / 10,      ... % sigma 1
       size(gaussImg, 1) / 10];         % sigma 2
% Define lower bound for parameter fitting
PLB = [0,   ...
       1,   ...
       1,   ...
       -1,  ...
       -90, ...
       0,   ...
       0];
% Define upper bound for parameter fitting
PUB = [2,             ...
       size(gaussImg, 1), ...
       size(gaussImg, 2), ...
       1,             ...
       90,            ...
       1e4,           ...
       1e4];    

% Fit to mean response to constrain subsequent fit center coords and ori
paramGssMu = lsqcurvefit('FuncWoNGaussian2DExpression', ...
                            PI, coord, dataMean, PLB, PUB);
                        
%% Function

function Z = FuncWoNGaussian2DExpression(P,Coordinates)
% 2-D Gaussian function WithOut Normalization
% Note: to create orientation, the whole coordinate is rotated.
%
% P = [Amp,X0,Y0,Base,Ort,SigmaMinor,SigmaMajor] --- parameters
% Coordinates.X --- X vector
% Coordinates.Y --- Y vector
% Z = Amp*exp(-(X-X0)^2/2/SigmaMajor^2-(Y-Y0)^2/2/SigmaMinor^2)+Base;

% Insert whatever code to define a 2-D Gaussian here
Z = [];

% Multiply the 2-D Gaussian by the expression map over space
Z = z.*Coordinates.expression;
end

