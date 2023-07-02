%% Define wanted sessions
% def sessions
sessions={'20230502','20230426'};

%% Load dataStruct
% def mainPath
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end

run([mainPath 'users/PK/colStimPipeline/GEVI/exptListGEVI.m'])
analysisSessID=[6]; %[7 4] %[6 5]  

%% Auto RunDA (add saTimeCourseROI)
open autoRunDA

%% Plot timecourse for ROI (requires ROITC.mat)
for sessionID=1:numel(analysisSessID)

    % define session IDs, session meta info
    currentSessID=analysisSessID(sessionID);
    dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto

    saveFlag=1;
    plotFlag=1;
    isSummary=1;

    analyzeTimecourse(dsCurrentSess)

end