function [adjustedSPD_uW, estimatedTotalPower_mW] = calculateSPD(behavioralData, imagingData, bitmapData, currentBlockStruct, imageNo, blockID, plotFlag)
    % Extract values from input struct
    pixelSize_mm = imagingData.pixelsizemm(blockID);
    pixelsON_bitmap = bitmapData.pixelsON(imageNo,blockID);
    pixelsON_DLP4710 = 5.8*10.4 / (pixelSize_mm)^2; % Total pixels in DMD rectangle
    timeON_ms = unique(behavioralData.TS.Header.Conditions.ProjTTLPulseOn);
    timeOFF_ms = unique(behavioralData.TS.Header.Conditions.ProjTTLPulseOff);
    cyclesON = unique(behavioralData.TS.Header.Conditions.ProjTTLnPulses);

    %% Calculate surface power density (accounting for pixels ON over illuminated DMD area)
    % Given power value across entire area
    estimatedTotalPower_mW=estimatePowerFromLED(currentBlockStruct, plotFlag); % total power (mW) over rectangular area recorded by the Thorlabs square sensor in mW
    
    % Calculate the fraction of the area illuminated
    fractionIlluminated = pixelsON_bitmap / pixelsON_DLP4710;
    
    % Adjust power delivered to the illuminated pixels
    powerDelivered_mW = estimatedTotalPower_mW * fractionIlluminated;
    
    % Calculate the area of one pixel
    pixelArea_mm2 = (pixelSize_mm)^2;
    
    % Calculate total illuminated area
    illuminatedArea_mm2 = pixelsON_bitmap * pixelArea_mm2;
    
    % Calculate Adjusted Surface Power Density (SPD)
    SPD = powerDelivered_mW / illuminatedArea_mm2; % mW/mm^2
      
    %% Calculate surface power density per second (accounting for pixels ON over illuminated DMD area)
    % Calculate Effective ON Time in seconds
    effectiveONTime_s = (timeON_ms * cyclesON) / 1000;
    effectiveOFFTime_s = (timeOFF_ms * cyclesON) / 1000;
    effectiveTotalTime_s = effectiveONTime_s + effectiveOFFTime_s;
    % Calculate Surface Power Density per Second (SPD_s)
    adjustedSPD_uW = SPD * 1000 * (effectiveONTime_s/effectiveTotalTime_s); % Normalize by effective total stimulation time (ON+OFF)
    
    % Display results
    fprintf(['Surface Power Density: %.2f ' char(181) 'W mm-2 (%.0f pixels)\n'], adjustedSPD_uW, pixelsON_bitmap);
end
