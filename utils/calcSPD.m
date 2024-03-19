function [SPD, SPD_s]=calcSPD(chamberID,LED,pixSize_mm,pixelsON,timeON_ms,cycles_ON,plotFlag)
% Calculates surface power density for given LED parameters and pixels in bitmap
% [SPD, SPD_s]=calcSPD(chamberID,LED,pixSize_mm,pixelsON,timeON_ms,cycles_ON,plotFlag)
% Inputs
pixSize_mm=behavioralData.Header.Conditions; % Pixel size in mm
%pixelsON= % Total pixels that are illuminated
estimatedPower_mW=estimatePowerFromLED(chamberID, LED, plotFlag) % Power (mW) recorded by the Thorlabs square sensor in mW
timeON_ms=behavioralData.Header.Conditions.ProjTTLPulseOn % LED ON time per cycle in ms
cycles_ON=behavioralData.Header.Conditions.ProjTTLnPulses % Total cycles of pulsing

% Calculate Pixel Area
pixelArea = pixSize_mm^2;

% Calculate Total Illuminated Area
illuminatedArea = pixelsON * pixelArea;

% Calculate Surface Power Density (SPD)
SPD = estimatedPower_mW / illuminatedArea;

% Calculate Effective ON Time in Seconds
effectiveTimeOn_s = (timeON_ms * cycles_ON) / 1000; % Effective ON time in seconds

% Calculate Surface Power Density per Second (SPD/s)
SPD_s = SPD / effectiveTimeOn_s;

% Display Results
fprintf('Illuminated Area: %f mm^2\n', illuminatedArea);
fprintf('Surface Power Density: %f mW/mm^2\n', SPD);
fprintf('Surface Power Density per Second: %f mW/mm^2/s\n', SPD_s);

end