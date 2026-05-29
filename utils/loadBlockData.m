function [currentBlockStruct,referenceBlockStruct,...
    blockData, behavioralData, imagingData, bitmapData, successFlag]=loadBlockData(datastruct, analysisBlockID, blockData, behavioralData, imagingData, bitmapData, blockID, pipelineMode, skipImaging)    
    disp('Loading behavioral, imaging and bitmap data...')
    %% Define filenames
    % Grab entries for blocks used
    [blockData, currentEntryID, currentEntry, referenceEntry, alignmentEntryDate] = getBlockInfo(datastruct, analysisBlockID, blockID, blockData);
    
    % Generate filenames for analyzed and reference block
    [currentBlockStruct,referenceBlockStruct]=grabFilenames(currentEntry, referenceEntry);

    % Get block paths
    [alignmentBlockPath] = getBlockPaths(currentEntry, alignmentEntryDate);
    
    %% Load behavioral, imaging and bitmap data
    % Load behav and imaging data
    [behavioralData, imagingData,successFlag]=loadBehavImagingData(currentBlockStruct, referenceBlockStruct, ...
        alignmentBlockPath, behavioralData, imagingData, blockID, pipelineMode, skipImaging);
    
    % Load bitmap params
    bitmapData = loadBitmapData(datastruct, currentBlockStruct, currentEntryID, bitmapData, blockID);
end