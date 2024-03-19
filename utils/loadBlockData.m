function [currentBlockStruct,referenceBlockStruct,...
    behavioralData, imagingData, bitmapData, successFlag]=loadBlockData(datastruct, analysisBlockID, behavioralData, imagingData, bitmapData, blockID)    
    disp('Loading behavioral, imaging and bitmap data...')
    %% Define filenames
    % Grab entries for blocks used
    [currentEntryID, currentEntry, referenceEntry, alignmentEntry] = getBlockInfo(datastruct,analysisBlockID,blockID);
    
    % Generate filenames for analyzed and reference block
    [currentBlockStruct,referenceBlockStruct]=grabFilenames(currentEntry, referenceEntry);

    % Get block paths
    [alignmentBlockPath] = getBlockPaths(alignmentEntry);
    
    %% Load behavioral, imaging and bitmap data
    % Load behav and imaging data
    [behavioralData, imagingData,successFlag]=loadBehavImagingData(currentBlockStruct, referenceBlockStruct, ...
        alignmentBlockPath, behavioralData, imagingData, blockID);
    
    % Load bitmap params
    bitmapData = loadBitmapData(datastruct, currentBlockStruct, currentEntryID, bitmapData, blockID);
    
    %
    
end