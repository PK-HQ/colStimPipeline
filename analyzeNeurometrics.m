%% Analyzes optostim data: 
% 1. Quality of VER (orientation map)
% 2. Bitmaps
% 3. Opto-evoked response

%% Load entries to dataStruct
%run('Y:\users\PK\colStimPipeline\exptListBiasing.m')
%datastructBias=datastruct;
run('Y:\users\PK\colStimPipeline\exptListBiasingFull.m')
datastructFull=datastruct;
%open('Y:\users\PK\colStimPipeline\exptListBiasing.m')
%open('Y:\users\PK\colStimPipeline\exptListBiasingFull.m')
%sessionIDs=[20221104, 20221109];
refSessID=38;
currSessID=39;
dataStructReference=datastructFull(refSessID);
dataStructCurrent=datastructFull(currSessID);

%% Assess orientation map quality (PCA explained variance, orientation map clarity)
%quality=compareOrtMaps(sessionIDs);

%% Assess bitmap quality (overlay of co-registered bitmaps)
%bitmaps=coregisterBitmaps(sessionIDs);

%% Assess opto effect (compare similarity of opto+visual with visual response)
assessOptoEffect2(dataStructReference,dataStructCurrent);
%assessOptoEffect(dataStructReference,dataStructCurrent);

