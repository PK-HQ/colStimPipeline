function [currentEntryID, currentEntry, referenceEntry, alignmentEntry] = getBlockInfo(datastruct,analysisBlockID,blockID)
    currentEntryID = analysisBlockID(blockID);
    currentEntry = datastruct(currentEntryID);
    referenceEntry = datastruct(currentEntry.referenceBlockNo);
    alignmentEntry = datastruct(currentEntry.alignmentBlockNo).date;
end