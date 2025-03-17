%% List of all reference, biasing and fixation blocks
% R-chamber: [12 14 15 17 19 20 21 23 24 26 29 31 32 34 37 40]; %17; %10 excluded because baseline different

%% Biasing sessions with entry no, points to session folder/run. Typically run 0 contains orientation map, run >=1 contains biasing experiment
% Baseline, Biasing (1x and 0.5x power, N=20), Fixation, PRF
entryNo=1; % O; % D+1 after cleaning
datastruct(entryNo).monkeyNo='32';
datastruct(entryNo).monkey='Pepper';
datastruct(entryNo).date='20250314';
datastruct(entryNo).run='0';
datastruct(entryNo).site='2';
datastruct(entryNo).modality='GCaMP';
datastruct(entryNo).alignmentBlockNo=1;
datastruct(entryNo).blueLED='L20OD10'; % lower to avoid hole in center
datastruct(entryNo).orangeLED='L30OD0';
datastruct(entryNo).greenImgID=[];
entryNo=2; % B;  Changed to gridsize 16
datastruct(entryNo).monkeyNo='32';
datastruct(entryNo).monkey='Pepper';
datastruct(entryNo).date='20250314';
datastruct(entryNo).run='3';
datastruct(entryNo).baselineTS='2'; %Run 2 corrupted
datastruct(entryNo).site='2';
datastruct(entryNo).modality='GCaMP';
datastruct(entryNo).alignmentBlockNo=1;
datastruct(entryNo).referenceBlockNo=1;
datastruct(entryNo).gridSize=16;%32;
datastruct(entryNo).gammaCorrFactor=[60 60];
datastruct(entryNo).sensitivity=.1*10^-1;
datastruct(entryNo).powercycle=4/50; %31-20 ON-OFF timing
datastruct(entryNo).blueLED='L20OD10'; % lower to avoid hole in center
datastruct(entryNo).orangeLED='L30OD0';
datastruct(entryNo).gaussianResponse='';
datastruct(entryNo).gaussianContourLevel=[];
datastruct(entryNo).gaussianContourLevelMax=[]; 
datastruct(entryNo).gaussianCond=[];
datastruct(entryNo).nColumnsWanted=[];
datastruct(entryNo).greenImgID=1;
