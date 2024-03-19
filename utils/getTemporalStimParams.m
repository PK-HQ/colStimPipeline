run([mainPath 'users/PK/colStimPipeline/exptListBiasingFull.m'])
for sessionID=1:numel(analysisSessID)%1:numel(analysisSessID)%1:numel(analysisSessID)%1:numel(analysisSessID)
    tic

    % define session IDs, session meta info
    currentSessID=analysisSessID(sessionID);
    dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
    filenameStructSession=generateFilenames(dsCurrentSess);
    load(filenameStructSession.TS);
    
    % Get temporal
    temporalTotal(sessionID)=unique(TS.Header.Conditions.ProjTTLPulseOn)+unique(TS.Header.Conditions.ProjTTLPulseOff);
    temporalDC(sessionID)=unique(TS.Header.Conditions.ProjTTLPulseOn)/temporalTotal(sessionID);
    temporalON(sessionID)=unique(TS.Header.Conditions.ProjTTLPulseOn)*unique(TS.Header.Conditions.ProjTTLnPulses);
end