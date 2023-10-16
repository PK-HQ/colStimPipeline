%% Debugging lost trials
clear all;

% Load biasing block
load('Y:\Chip\Chip20230816\run1\M28D20230816R1TS.mat')

% Grab outcome codes and counts for trials within block
trialOutcomes=vertcat(TS.Trial.Outcome);
uniqueOutcomeIDs=double(unique(trialOutcomes));
uniqueOutcomeIDs=uniqueOutcomeIDs(abs(uniqueOutcomeIDs)>=30);
cont(:,1)=uniqueOutcomeIDs;
for id=1:length(uniqueOutcomeIDs)
    count=numel(find(trialOutcomes==uniqueOutcomeIDs(id)));
    cont(id,2)=count;
end
% calculate percentage of error trials in block
cont(:,3)=round(cont(:,2)*100./numel(trialOutcomes),0);

% Show error codes and error count
%{
   NOT_ACQUIRE_TARG: -52
      NOT_HOLD_TARG: -53
        NOT_HOLD_FA: -54
NOT_HOLD_WRONG_TARG: -55
NOT_HOLD_IRREL_TARG: -56
      SACCADE_EARLY: -57
%}
t = array2table(cont,'VariableNames',{'Error code' 'Count' '% of trials'});
disp(t)
disp(TS.Header.DEF.OUTCOME)
