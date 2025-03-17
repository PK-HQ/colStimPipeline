function [str,isMerged]=checkMergedBaseline(datastruct,analysisBlockID,blockID)
isMerged=isempty(datastruct(analysisBlockID(blockID)).baselineTS);
switch isMerged
    case 0
        str='split';
    case 1
        str='combined';
end