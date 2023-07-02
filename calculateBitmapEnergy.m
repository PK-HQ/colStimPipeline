function bitmapEnergy=calculateBitmapEnergy(TS, orangeLightPD, pixelDensity)
% orangeLightPDpower density measured with power sensor
spatialPower=orangeLightPD*pixelDensity/100; %in mW
pulseDuration=TS.Header.ConditionParams.Projector.Duration_On__ms*TS.Header.ConditionParams.Projector.Number_of_Pulses / 1000; %in seconds
bitmapEnergy = spatialPower * pulseDuration; %mW/s
end