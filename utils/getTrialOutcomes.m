function [outcomeMat, outcomeTbl]=getTrialOutcomes(TSstruct)
%% Get outcomes per trial (correct/error etc)
%here are outcome defs
outcomeDefinitions=TSstruct.TS.Header.DEF.OUTCOME;
successCode=10;
failedCode=11;

%get outcome of each trial
outcomeTrials=extractfield(TSstruct.TS.Trial, 'Outcome');

%get correct and error as binary array
correctTrialwise=double(outcomeTrials==successCode);
errorTrialwise=double(outcomeTrials==failedCode);
errorTrialwise(errorTrialwise==1)=-1;

%code correct=1, error=0, others=NaN
allTrialwise=correctTrialwise+errorTrialwise;
allTrialwise(allTrialwise==0)=NaN;
allTrialwise(allTrialwise==-1)=0;

% init
nTrials=size(outcomeTrials,2);
gaborOrt=NaN(1,nTrials);
projImg=NaN(1,nTrials);
congruentVisOptostim=NaN(1,nTrials);

%% Get stimulus presented per trial
for trialNo=1:size(outcomeTrials,2)
    trialGraphics=TSstruct.TS.Trial(trialNo).Graphics;
    if size(trialGraphics,2)>4 && contains(trialGraphics{1,4}.Filename,'.bmp')
        %ort and contrast, ort = orientation * contrast (H=-1,V=+1)
        gaborOrtTrial=trialGraphics{1,4}.Filename; %Oxxxx
        stimConTrial=trialGraphics{1,4}.Gain * 100;
        gaborOrtTrials(trialNo)=trialGraphics{1,4}.Filename;
        if contains(gaborOrtTrial,'O00000')
            gaborOrt(trialNo)=-1 * abs(stimConTrial); %H
        elseif  contains(gaborOrtTrial,'O09000')
            gaborOrt(trialNo)=1 * abs(stimConTrial); %V
        end
        %recode projector image as -1, 0, 1
        projImgTrial=TSstruct.TS.Trial((trialNo)).Graphics{1,3}.Filename; %last \Dot OR \Oxxx.bmp
        if contains(projImgTrial,'Dot')
            projImg(trialNo)=0; %baseline, black bmp stim
        elseif contains(projImgTrial,'O00000')
            projImg(trialNo)=-1; %H
        elseif  contains(projImgTrial,'O09000')
            projImg(trialNo)=1; %V
        else
            %projImg(trialNo)=0; %baseline, nostim (for 20221028)
        end
        
        if sign(gaborOrt(trialNo)) == sign(projImg(trialNo))
            congruentVisOptostim(trialNo)=1;
        elseif sign(gaborOrt(trialNo)) == -1*sign(projImg(trialNo))
            congruentVisOptostim(trialNo)=-1;
        elseif sign(projImg(trialNo))==0
            congruentVisOptostim(trialNo)=0;
        end
    else
        gaborOrt(trialNo)=NaN;
        projImg(trialNo)=NaN;
        congruentVisOptostim(trialNo)=NaN;
    end
end




%% Make into a mat and table
%keep only valid trials
validTrials=find(~isnan(allTrialwise) & ~isnan(gaborOrt)); %is valid trial with correct/error & non-blank
validTrialNo=1:numel(validTrials);

%filter
gaborOrt=gaborOrt(validTrials); %NaN = blanks
projImg=projImg(validTrials); %NaN = blanks
congruentVisOptostim=congruentVisOptostim(validTrials); %NaN = blanks
allTrialwise=allTrialwise(validTrials);

%make tbl and mat
headers={'ValidTrialCount','Visualstim','Optostim','CongruentVO','Correct'};
outcomeTbl=table(validTrialNo',gaborOrt',projImg',congruentVisOptostim',allTrialwise','VariableNames',headers);
outcomeMat=[validTrialNo',gaborOrt',projImg',congruentVisOptostim',allTrialwise'];

%print first 10 rows
disp(outcomeTbl(:,:))
end