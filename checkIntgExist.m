clear all;clc

%% From mainColstimPipeline
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end
saveFlag=0;
run([mainPath 'users/PK/colStimPipeline/exptListBiasingFull.m'])

analysisSessID=fliplr([93 92 90 89 88 87 84 83 81 79 78 74 73 72 67 66 65 64 62 59 51 50])

for sessionID=1:numel(analysisSessID)
    currentSessID=analysisSessID(sessionID);
    
    dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
    referenceSessID=dsCurrentSess.referenceSessionEntryNo;
    dsReferenceSess=datastruct(referenceSessID); %fixed, to get ort map projected onto

    filenameStructCurrent=generateFilenames(dsCurrentSess);
    filenameStructReference=generateFilenames(dsReferenceSess);

    %% Stage 1. Load integrated cortical response, TS, gaussian ROI mask
    if isfile(filenameStructCurrent.Intg) % Load imaging data if it exists 
    else
        disp(sprintf('Missing Intg.: #%.0f',currentSessID))
    end
end