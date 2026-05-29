function [meanPowerDensityWithinROI_mWmm2, totalPowerToOnPixelsWithinROI_mW, ...
          projectorPowerDensity_mWmm2, temporalDutyCycle] = calculateSPD( ...
          behavioralData, imagingData, bitmapData, currentBlockStruct, imageNo, blockID, plotFlag)

    % Full-white projector optical power density.
    % Units: mW/mm^2
    projectorPowerDensity_mWmm2 = estimatePowerFromLED(currentBlockStruct, bitmapData, plotFlag);

    % Pull spatial terms already computed in convertForProjectorGPT.
    % spatialDutyCycleWithinROI is unitless.
    % areaPixelsONWithinROI is mm^2.
    spatialDutyCycleWithinROI = bitmapData.spatialDutyCycleWithinROI(imageNo, blockID);
    pixelsOnAreaWithinROI_mm2 = bitmapData.areaPixelsONWithinROI(imageNo, blockID);

    % Get pulse timing.
    if ~isempty(behavioralData)
        timeON_ms = unique(behavioralData.optoTS(blockID).Header.Conditions.ProjTTLPulseOn);
        timeOFF_ms = unique(behavioralData.optoTS(blockID).Header.Conditions.ProjTTLPulseOff);
    else
        timeON_ms = bitmapData.ProjTTLPulseOn(blockID);
        timeOFF_ms = bitmapData.ProjTTLPulseOff(blockID);
    end

    % Temporal duty cycle within the pulse cycle.
    % Example: 25 ms ON, 25 ms OFF -> 0.5.
    temporalDutyCycle = timeON_ms ./ (timeON_ms + timeOFF_ms);

    % Boss-style metric:
    % Mean optical power density over the final ROI and stimulation window.
    % Units: mW/mm^2
    meanPowerDensityWithinROI_mWmm2 = ...
        projectorPowerDensity_mWmm2 .* spatialDutyCycleWithinROI .* temporalDutyCycle;

    % Old-style absolute-power audit metric:
    % Total average optical power delivered to ON pixels inside the final ROI.
    % Units: mW
    totalPowerToOnPixelsWithinROI_mW = ...
        projectorPowerDensity_mWmm2 .* pixelsOnAreaWithinROI_mm2 .* temporalDutyCycle;

    % Display results.
    fprintf('Projector power density: %.4f mW/mm^2\n', projectorPowerDensity_mWmm2);
    fprintf('Spatial duty cycle within ROI: %.4f\n', spatialDutyCycleWithinROI);
    fprintf('Temporal duty cycle: %.4f\n', temporalDutyCycle);
    fprintf('Mean power density within ROI: %.4f mW/mm^2\n', meanPowerDensityWithinROI_mWmm2);
    fprintf('Total power to ON pixels within ROI: %.4f mW\n', totalPowerToOnPixelsWithinROI_mW);
end