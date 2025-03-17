function [alignmentBlockPath] = getBlockPaths(currentEntry, alignmentBlock)
    computerName = getenv('COMPUTERNAME');
    hostname = getenv('HOSTNAME');
    alignmentBlockPath = '';
    %currentBlockPath = '';
    monkeyName=currentEntry.monkey;
    if ispc
        switch computerName
            case {'LA-CPSD077020WD', 'LA-CPSA07019WD', 'CVIS-A64882', 'PSYC-A77304'}
                alignmentBlockPath = ['Y:/' monkeyName '/' monkeyName alignmentBlock '/'];
                if ismember(computerName, {'CVIS-A64882', 'PSYC-A77304'})
                    %currentBlockPath = ['Y:/' 'Chip/Chip' currentBlock.date '/'];
                end
        end
    elseif contains(hostname, 'psy.utexas.edu')
        alignmentBlockPath = ['/eslab/data/'  monkeyName '/' monkeyName alignmentBlock '/'];
        %currentBlockPath = ['/eslab/data/' 'Chip/Chip' currentBlock.date '/'];
    end
end