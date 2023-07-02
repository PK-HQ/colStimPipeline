function [trials,images,condIDs]=getUsableTrials(TS,ImagingData)
%% [README] Extracts completed trials from TS
% Input:
% TS - TS file with metadata for each trial
% ImagingData - trial-wise imaging data file

% Output: 
% trials - struct that contains all completed trials, correct trials and error trials
% images - struct that contains cortical images for the above, and also averaged across all completed
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

% get integrated images of usable trials
if ndims(ImagingData)==3
    images.All=[];
    images.Average=ImagingData(:,:,usableTrialIdx);
    images.Correct=[];
    images.Error=[];    
elseif ndims(ImagingData)==4
    images.All=ImagingData(:,:,:,usableTrialIdx);
    images.Average=ImagingData(:,:,usableTrialIdx);
    images.Correct=ImagingData(:,:,:,usableCorrTrialIdx);
    images.Error=ImagingData(:,:,:,usableErrTrialIdx);    
end



nTrials=numel(trials.All);
nImages=size(images.Average,3);

if nTrials == nImages
    % nothing
else
   error('Error: No. of trials and images do not match.') 
end

condIDs.blankConds = find(cellfun(@isempty,TS.Header.Conditions.ProjImg));
condIDs.baselineConds = findMatch(TS.Header.Conditions.ProjImg,'Dot');
condIDs.opto0Conds = findMatch(TS.Header.Conditions.ProjImg,'O000');
condIDs.opto90Conds = findMatch(TS.Header.Conditions.ProjImg,'O090');
end