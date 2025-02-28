function [energy, powerdensity, timeONPercent] = calculateSPD(behavioralData, imagingData, bitmapData, currentBlockStruct, imageNo, blockID, plotFlag)
method=2;
switch method
    case {1}
    % Extract values from input struct
    pixelSize_mm = imagingData.pixelsizemm(blockID);
    pixelsON_bitmap = bitmapData.pixelsON(imageNo,blockID);
    pixelsON_DLP4710 = 5.8*10.4 / (pixelSize_mm)^2; % Total pixels in DMD rectangle
    timeON_ms = unique(behavioralData.optoTS(blockID).Header.Conditions.ProjTTLPulseOn);
    timeOFF_ms = unique(behavioralData.optoTS(blockID).Header.Conditions.ProjTTLPulseOff);
    timeONPercent = timeON_ms/(timeON_ms+timeOFF_ms);
    cyclesON = unique(behavioralData.optoTS(blockID).Header.Conditions.ProjTTLnPulses);

    %% Calculate surface power density (accounting for pixels ON over illuminated DMD area)
    % Given power value across entire area
    powerdensity=estimatePowerFromLED(currentBlockStruct, plotFlag); % total power (mW) over rectangular area recorded by the Thorlabs square sensor in mW
    % Calculate the fraction of the area illuminated
    fractionIlluminated = pixelsON_bitmap / pixelsON_DLP4710;
    % Adjust power delivered to the illuminated pixels
    powerDelivered_mW = powerdensity * fractionIlluminated;
    
    % Calculate the area of one pixel
    pixelArea_mm2 = (pixelSize_mm)^2;
    % Calculate total illuminated area
    illuminatedArea_mm2 = pixelsON_bitmap * pixelArea_mm2;
    
    % Calculate Adjusted Surface Power Density (SPD)
    energy = powerDelivered_mW / illuminatedArea_mm2; % mW/mm^2
      
    %% Calculate surface power density per second (accounting for pixels ON over illuminated DMD area)
    % Calculate Effective ON Time in seconds
    onTime_s = (timeON_ms * cyclesON) / 1000;
    offTime_s = (timeOFF_ms * cyclesON) / 1000;
    totalTime_s = onTime_s + offTime_s;
    % Calculate Surface Power Density per Second (SPD_s)
    adjustedSPD_uW = energy * 1000 * (onTime_s/totalTime_s); % Normalize by effective total stimulation time (ON+OFF)
    
    case {2}
        % Extract values from input struct
        powerdensity=estimatePowerFromLED(currentBlockStruct, plotFlag); %mW/mm2, total power (mW) over rectangular area recorded by the Thorlabs square sensor in mW
        pixelsOnArea= bitmapData.areaPixelsON(imageNo,blockID); %mm2
        power = powerdensity * pixelsOnArea; %mW
        
        timeON_ms = unique(behavioralData.optoTS(blockID).Header.Conditions.ProjTTLPulseOn);
        timeOFF_ms = unique(behavioralData.optoTS(blockID).Header.Conditions.ProjTTLPulseOff);
        timeONPercent = timeON_ms/(timeON_ms+timeOFF_ms);
        
        energy = power * timeONPercent;
end
    % Display results
    fprintf(['Energy delivered: %.2f mW over 300ms\n'], energy);
end
