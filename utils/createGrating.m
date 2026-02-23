%% Create square sine-grating BMP files (no Gaussian fall-off)
% YC at ES lab – revised 22-Aug-2025

%% Clear workspace
clearvars; close all;

%% Parameters you supplied -----------------------------------------------
RatioPxlPerDeg = 50;                       % pixels / deg

refPos   = [2 -3.5];
refSize  = 2;
refSF    = 2;
currentPos = [-1.1 -2.2];

[newSize,newSf] = scaleStimByEcc(refPos,refSize,refSF,currentPos,'mean');

GaborSizes = 4; %round(newSize,1);             % deg
GaborSFs   = 3.8; %round(newSf ,1);              % cpd
GaborOrts  = [0:15:180];                       % deg
% -------------------------------------------------------------------------

%% Where to save
PathName = uigetdir('.','Select directory to save BMP files');
if isequal(PathName,0); disp('No directory selected.'); return; end

%% Generate stimuli
for GaborSize = GaborSizes                               % deg (square)
    
    SizePix = ceil(GaborSize*RatioPxlPerDeg/2)*2 + 1;    % odd # pixels
    [X,Y]   = meshgrid( -(SizePix-1)/2 : (SizePix-1)/2 );% pixel coords
    
    % Blank (mid-gray)
    imwrite(0.5*ones(SizePix), ...
        fullfile(PathName, sprintf('BlankZ%04g.bmp',GaborSize*100)), 'bmp');
    
    % Gratings
    for GaborOrt = GaborOrts                              % deg
        theta = deg2rad(GaborOrt);                        % rad
        Xrot  =  X*cos(theta) + Y*sin(theta);             % axis orth. to bars
        
        for GaborPhs = 0                                  % deg (change as needed)
            phaseRad = deg2rad(GaborPhs);                 % rad
            
            for GaborSF = GaborSFs                        % cpd
                sf_pix = GaborSF / RatioPxlPerDeg;        % cycles / pixel
                
                % Pure sinusoidal grating, range [0,1]
                A = 0.5 + 0.5 * cos( 2*pi*sf_pix*Xrot + phaseRad );
                
                % Save
                imwrite(A, fullfile(PathName, ...
                    sprintf('SGZ%04gS%04dO%05dP%05d.bmp', ...
                    GaborSize*100, GaborSF*100, GaborOrt*100, GaborPhs*100)), ...
                    'bmp');
            end
        end
    end
end

%% ---------------- helper (unchanged) ------------------------------------
function [newSize,newSf] = scaleStimByEcc(prevPos,prevSize,prevSf,newPos,cycleMode)
if nargin < 5, cycleMode = 'mean'; end
eccPrev = hypot(prevPos(:,1),prevPos(:,2));
eccNew  = hypot(newPos(1), newPos(2));
sizePerEcc = prevSize./eccPrev;
newSize    = max(sizePerEcc) * eccNew;
cyclesPrev = prevSize .* prevSf;
switch lower(cycleMode)
    case 'min',  tgtCycles = min(cyclesPrev);
    case 'max',  tgtCycles = max(cyclesPrev);
    otherwise,   tgtCycles = mean(cyclesPrev);
end
newSf = tgtCycles / newSize;
end
