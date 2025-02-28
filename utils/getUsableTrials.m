function [trials,images,condIDs]=getUsableTrials(TS, ImagingData)
%% [README] Extracts completed trials from TS
% Input:
% TS - TS file with metadata for each trial
% ImagingData - trial-wise imaging data file

% Output: 
% trials - struct that contains all completed trials, correct trials and error trials
% images - struct that contains cortical images for the above, and also averaged across all completed
%           - images.trials  = 512 x 512 x 32 conditions x 10 trials
%           - images.average = 512 x 512 x 32 conditions x 1 average across trials
% condIDs - condition IDs for blank, baseline, opto0, opto90 conds (e.g. [1 3 5])

% Define relevant allegro code IDs
allegroCodingScheme.successCode=10;
allegroCodingScheme.failedCode=11;
allegroCodingScheme.successSpecialCode=-10;
allegroCodingScheme.failedSpecialCode=-11;

% Get fields that contain outcome ID (Outcome), condition ID (CurrCond) and trials with imaging collected (FlagOIBLK)
trialOutcomeID=double(extractfield(TS.Trial, 'Outcome'));
trialCondID=double(extractfield(TS.Trial, 'CurrCond'));
trialOIID=double(extractfield(TS.Trial, 'FlagOIBLK'));
trialOIID(trialOIID==0)=NaN;

% Get trialID of each outcome type (with imaging, correct and error)
trialOutcomeID=trialOutcomeID.*trialOIID;% get trials with saved image, those without are set to nan (FlagOIBLK==1)
completedTrialIdx=find(trialOutcomeID == allegroCodingScheme.successCode | trialOutcomeID == allegroCodingScheme.failedCode |...
    trialOutcomeID == allegroCodingScheme.successSpecialCode | trialOutcomeID == allegroCodingScheme.failedSpecialCode)'; % get completed trials (Outcome==10,11,-10,-11)
  
% get index of all correct and error trials
usableTrialIdxTmp=find(trialOutcomeID == allegroCodingScheme.successCode | trialOutcomeID == allegroCodingScheme.failedCode); % filter for Outcome==10,11 and indexes it to be within 1:nTotalTrials, e.g. 1:363
% get index of all correct trials
usableCorrTrialIdxTmp=find(trialOutcomeID == allegroCodingScheme.successCode); % indexes it to be within 1:nTotalTrials, e.g. 1:363
% get index of all error trials
usableErrTrialIdxTmp=find(trialOutcomeID == allegroCodingScheme.failedCode); % indexes it to be within 1:nTotalTrials, e.g. 1:363

% get usable trials, with trial index based out of completed trials
[usableTrialIdx,~]=find(completedTrialIdx==usableTrialIdxTmp); % re-indexes index of 1:nTotalTrials to be within 1:nCompletedTrials, e.g. 1:266
[usableCorrTrialIdx,~]=find(completedTrialIdx==usableCorrTrialIdxTmp); % re-indexes index of 1:nTotalTrials to be within 1:nCompletedTrials, e.g. 1:266
[usableErrTrialIdx,~]=find(completedTrialIdx==usableErrTrialIdxTmp); % re-indexes index of 1:nTotalTrials to be within 1:nCompletedTrials, e.g. 1:266

% get condition of each trial
trials.All=trialCondID(usableTrialIdxTmp)'; % grabs trials using index of 1:nTotalTrials
trials.Correct=trialCondID(usableCorrTrialIdxTmp)'; % grabs trials using index of 1:nTotalTrials
trials.Error=trialCondID(usableErrTrialIdxTmp)'; % grabs trials using index of 1:nTotalTrials

%% Split imaging data into correct, error and combined
if ndims(ImagingData)==3
    
    %% Sort condIDs into visual x opto conditions
    % Find conditions that are blank, baseline (visual only), opto0, opto90
    condIDs.blankConds = find(cellfun(@isempty,TS.Header.Conditions.ProjImg)); nBlanks=numel(condIDs.blankConds);
    condIDs.baselineConds = findMatch(TS.Header.Conditions.ProjImg,'Dot');
    condIDs.opto0Conds = findMatch(TS.Header.Conditions.ProjImg,'O000');
    condIDs.opto90Conds = findMatch(TS.Header.Conditions.ProjImg,'O090');
    % find conditions where visual gabor is 0 and 90
    gabor0=find(TS.Header.Conditions.GaborOrt==0 & TS.Header.Conditions.TypeCond == 3);
    gabor90=find(TS.Header.Conditions.GaborOrt==90 & TS.Header.Conditions.TypeCond == 3);
    % Sort the above into the following visual and opto pairings
    condIDs.V0=condIDs.baselineConds(ismember(condIDs.baselineConds, gabor0)); %condIDs.baselineConds(1:numel(condIDs.baselineConds)/2);%[1:3:nOptoConds/2]+nblank; 
    condIDs.V90=condIDs.baselineConds(ismember(condIDs.baselineConds, gabor90)); %condIDs.baselineConds(numel(condIDs.baselineConds)/2 + 1 : numel(condIDs.baselineConds));%[nOptoConds/2+1:3:nOptoConds]+nblank;
    condIDs.V0O0=condIDs.opto0Conds(ismember(condIDs.opto0Conds, gabor0)); %condIDs.opto0Conds(1:numel(condIDs.opto0Conds)/2);
    condIDs.V90O0=condIDs.opto0Conds(ismember(condIDs.opto0Conds, gabor90));
    condIDs.V0O90=condIDs.opto90Conds(ismember(condIDs.opto90Conds, gabor0));
    condIDs.V90O90=condIDs.opto90Conds(ismember(condIDs.opto90Conds, gabor90));

    %% Subtract blank, remove blank
    [blankTrialIdx,~]=find(trials.All==condIDs.blankConds);
    blankTrialsAvg=mean(ImagingData(:,:,blankTrialIdx),3);
    ImagingDataBlankSubt=ImagingData-blankTrialsAvg;
    
    %% Average the imaging data by condition
    exptConds=1:max(trials.All);
    exptConds(condIDs.blankConds)=[];
    for condNo=1:numel(exptConds)
        condID=exptConds(condNo);
        % get trial indices
        [condTrialIdx,~]=find(trials.All==condID);
        allTrialIdx=intersect([usableCorrTrialIdx;usableErrTrialIdx],condTrialIdx); % EDIT added 10/29 to combine correct and error 
        correctTrialIdx=intersect(usableCorrTrialIdx,condTrialIdx);
        errorTrialIdx=intersect(usableErrTrialIdx,condTrialIdx);

        % store trial images by condition x correct/error
        images.trials(:,:,condNo,1:numel(allTrialIdx))=ImagingDataBlankSubt(:,:,allTrialIdx);
        images.trialsCorrect(:,:,condNo,1:numel(correctTrialIdx))=ImagingDataBlankSubt(:,:,correctTrialIdx);
        images.trialsError(:,:,condNo,1:numel(errorTrialIdx))=ImagingDataBlankSubt(:,:,errorTrialIdx);
        
        % store average images by condition x correct/error
        images.average(:,:,condNo)=nanmean(images.trials(:,:,condNo,:),4);
        images.averageCorrect(:,:,condNo)=nanmean(images.trialsCorrect(:,:,condNo,:),4);
        images.averageError(:,:,condNo)=nanmean(images.trialsError(:,:,condNo,:),4);
    end
   
   % Remove all blank conds from idx
   fields=fieldnames(condIDs);
   for fn=1:numel(fields)
       name=fields{fn};
      % reassign
       condIDs.(name)=condIDs.(name)-nBlanks; %=condIDs.(fn{1})
   end
   
elseif ndims(ImagingData)==4
    %{
    images.All=ImagingDataBlankSubt(:,:,:,usableTrialIdx);
    images.Average=ImagingDataBlankSubt(:,:,usableTrialIdx);
    images.Correct=ImagingDataBlankSubt(:,:,:,usableCorrTrialIdx);
    images.Error=ImagingDataBlankSubt(:,:,:,usableErrTrialIdx);    
    %}
end
end