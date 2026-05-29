function [blockData, currentEntryID, currentEntry, referenceEntry, alignmentEntryDate] = getBlockInfo(datastruct,analysisBlockID ,blockID, blockData)
    currentEntryID = analysisBlockID(blockID);
    currentEntry = datastruct(currentEntryID);
    referenceEntry = datastruct(currentEntry.referenceBlockNo);
    alignmentEntry = datastruct(currentEntry.alignmentBlockNo);

    alignmentEntryDate = alignmentEntry.date;

    %% Store block data
    %file for alignment
    blockData.alignment_date{blockID}=alignmentEntry.date;
    blockData.alignment_run{blockID}=alignmentEntry.run; 
    %file for reference ort map
    blockData.ortmap_date{blockID}=referenceEntry.date;
    blockData.ortmap_run{blockID}=referenceEntry.run; 
    %file for gaussian ort map
    blockData.gaussfootprint_run{blockID}=currentEntry.date;
    blockData.gaussfootprint_run{blockID}=currentEntry.gaussianResponse; 
     % file for baseline block
    blockData.baseline_date{blockID}=currentEntry.date;
    blockData.baseline_run{blockID}=currentEntry.baselineTS;
end