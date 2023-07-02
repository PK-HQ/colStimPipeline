function quality=compareOrientationMaps(sessionIDs)
%% ----- Get params ----- 
nSessions=numel(sessionIDs);

%% ----- Load data across sessions -----
for sessionNo=1:nSessions
    % generate filenames
    filename=genFilenames();
    
    % load imaging files into container
    multiSessData(sessionNo).ortmap=load(); % imaging files for passive viewing of orientations
    
    % load green images into container
    multiSessData(sessionNo).greenimg=load(); % imaging files for passive viewing of orientations

end


%% ----- Calculate single session orientation map quality statistics ----- 
for sessionNo=1:nSessions
    sessData=multiSessData(sessionNo);
    quality(sessionNo,:)=calcOrtMapQC(sessData);
end

function calcOrtMapQC()
    % get orientation map with amplitude
    % get tuning curve sigma
    % get PCA explained var per component
    % get overlap between 0 and 90
end

%% ----- Calculate cross-session orientation map quality statistics ----- 

end