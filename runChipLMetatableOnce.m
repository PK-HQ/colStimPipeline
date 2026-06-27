function runChipLMetatableOnce()
% Compatibility wrapper — regenerate Chip L metatable.
% Delegates all logic to runMetatableForTarget.
%
% Usage:
%   cd('Y:\users\PK\colStimPipeline')
%   runChipLMetatableOnce

restoredefaultpath;
addpath(genpath('Y:/users/PK/colStimPipeline'));

result = runMetatableForTarget('Y:/', 'Chip', 'L');
if result.success
    fprintf('Chip L metatable complete: %s\n', result.matOutputPath);
else
    error('runChipLMetatableOnce:Failed', ...
        'runMetatableForTarget returned success=false for Chip L');
end
end
