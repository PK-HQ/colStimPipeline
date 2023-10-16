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
%get trials with OI BLK
trialOutcomeID=trialOutcomeID.*trialOIID;
completedTrialIdx=find(trialOutcomeID == allegroCodingScheme.successCode | trialOutcomeID == allegroCodingScheme.failedCode |...
    trialOutcomeID == allegroCodingScheme.successSpecialCode | trialOutcomeID == allegroCodingScheme.failedSpecialCode)'; 
  
% get index of all correct and error trials
usableTrialIdxTmp=find(trialOutcomeID == allegroCodingScheme.successCode | trialOutcomeID == allegroCodingScheme.failedCode); % indexes it to be within 1:nTotalTrials, e.g. 1:363
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
    % Grouped by optostim condition
    condIDs.blankConds = find(cellfun(@isempty,TS.Header.Conditions.ProjImg));
    condIDs.baselineConds = findMatch(TS.Header.Conditions.ProjImg,'Dot');
    condIDs.opto0Conds = findMatch(TS.Header.Conditions.ProjImg,'O000');
    condIDs.opto90Conds = findMatch(TS.Header.Conditions.ProjImg,'O090');
    % Grouped by visual x optostim condition
    condIDs.V0=condIDs.baselineConds(1:numel(condIDs.baselineConds)/2);%[1:3:nOptoConds/2]+nblank; 
    condIDs.V90=condIDs.baselineConds(numel(condIDs.baselineConds)/2 + 1 : numel(condIDs.baselineConds));%[nOptoConds/2+1:3:nOptoConds]+nblank;
    condIDs.V0O0=condIDs.opto0Conds(1:numel(condIDs.opto0Conds)/2);
    condIDs.V90O0=condIDs.opto0Conds(numel(condIDs.opto0Conds)/2 + 1 : numel(condIDs.opto0Conds));
    condIDs.V0O90=condIDs.opto90Conds(1:numel(condIDs.opto90Conds)/2);
    condIDs.V90O90=condIDs.opto90Conds(numel(condIDs.opto90Conds)/2 + 1 : numel(condIDs.opto90Conds));

    %% Subtract blank
    [blankTrialIdx,~]=find(trials.All==condIDs.blankConds);
    blankTrialsAvg=mean(ImagingData(:,:,blankTrialIdx),3);
    ImagingDataBlankSubt=ImagingData-blankTrialsAvg;
    
    images.trials=nan(size(ImagingDataBlankSubt,1),size(ImagingDataBlankSubt,2),max(trials.All),20);
    images.trialsCorrect=nan(size(ImagingDataBlankSubt,1),size(ImagingDataBlankSubt,2),max(trials.All),20);
    images.trialsError=nan(size(ImagingDataBlankSubt,1),size(ImagingDataBlankSubt,2),max(trials.All),20);

    images.average=nan(size(ImagingDataBlankSubt,1),size(ImagingDataBlankSubt,2),max(trials.All));
    images.averageCorrect=nan(size(ImagingDataBlankSubt,1),size(ImagingDataBlankSubt,2),max(trials.All));
    images.averageError=nan(size(ImagingDataBlankSubt,1),size(ImagingDataBlankSubt,2),max(trials.All));
    
    %% Average the imaging data by condition
    for condID=1:max(trials.All)
        % get trial indices
        [condTrialIdx,~]=find(trials.All==condID);
        correctTrialIdx=intersect(usableCorrTrialIdx,condTrialIdx);
        errorTrialIdx=intersect(usableErrTrialIdx,condTrialIdx);

        % store trial images by condition x correct/error
        images.trials(:,:,condID,1:numel(condTrialIdx))=ImagingDataBlankSubt(:,:,condTrialIdx);
        images.trialsCorrect(:,:,condID,1:numel(correctTrialIdx))=ImagingDataBlankSubt(:,:,correctTrialIdx);
        images.trialsError(:,:,condID,1:numel(errorTrialIdx))=ImagingDataBlankSubt(:,:,errorTrialIdx);
        
        % store average images by condition x correct/error
        images.average(:,:,condID)=nanmean(images.trials(:,:,condID,:),4);
        images.averageCorrect(:,:,condID)=nanmean(images.trialsCorrect(:,:,condID,:),4);
        images.averageError(:,:,condID)=nanmean(images.trialsError(:,:,condID,:),4);
    end
    %images.Average=ImagingDataBlankSubt(:,:,usableTrialIdx);
elseif ndims(ImagingData)==4
    %{
    images.All=ImagingDataBlankSubt(:,:,:,usableTrialIdx);
    images.Average=ImagingDataBlankSubt(:,:,usableTrialIdx);
    images.Correct=ImagingDataBlankSubt(:,:,:,usableCorrTrialIdx);
    images.Error=ImagingDataBlankSubt(:,:,:,usableErrTrialIdx);    
    %}
end
end