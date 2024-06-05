
[mainPath, datastruct]=setupEnv('users/PK/colStimPipeline/exptListBiasingFull.m');
nColumnsWanted=[]; chamberWanted='R';
analysisBlockID = organizeBlocks(datastruct, chamberWanted, nColumnsWanted);
analysisBlockID_R=analysisBlockID;
tblR=data_table(analysisBlockID_R,:);
table2excel(tblR, [mainPath 'users\PK\colStimPipeline\docs\optostimTableR'])

nColumnsWanted=[]; chamberWanted='L';
analysisBlockID = organizeBlocks(datastruct, chamberWanted, nColumnsWanted);
analysisBlockID_L=analysisBlockID;
tblL=data_table(analysisBlockID_L,:);
table2excel(tblL, [mainPath 'users\PK\colStimPipeline\docs\optostimTableL'])
