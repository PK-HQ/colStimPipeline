function  [V0, V90, V0O0, V90O0, V0O90, V90O90]=getCondIDs(condIDs)
% Define condition indices
V0=condIDs.baselineConds(1:numel(condIDs.baselineConds)/2);%[1:3:nOptoConds/2]+nblank; 
V90=condIDs.baselineConds(numel(condIDs.baselineConds)/2 + 1 : numel(condIDs.baselineConds));%[nOptoConds/2+1:3:nOptoConds]+nblank;
V0O0=condIDs.opto0Conds(1:numel(condIDs.opto0Conds)/2);
V90O0=condIDs.opto0Conds(numel(condIDs.opto0Conds)/2 + 1 : numel(condIDs.opto0Conds));
V0O90=condIDs.opto90Conds(1:numel(condIDs.opto90Conds)/2);
V90O90=condIDs.opto90Conds(numel(condIDs.opto90Conds)/2 + 1 : numel(condIDs.opto90Conds));
end