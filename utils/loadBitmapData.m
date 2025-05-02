function bitmapData = loadBitmapData(datastruct, currentBlockStruct, currentBlockID, bitmapData, blockID)
    bitmapData.gridSize(:,blockID) = datastruct(currentBlockID).gridSize;
    if length(datastruct(currentBlockID).gammaCorrFactor)==1
        datastruct(currentBlockID).gammaCorrFactor=repmat(datastruct(currentBlockID).gammaCorrFactor,1,2);
    end
    
    bitmapData.gammaCorrFactor(:,:,blockID) = datastruct(currentBlockID).gammaCorrFactor;
    
    newThresholdFn=~isnan(datastruct(currentBlockID).sensitivity);
    switch newThresholdFn
        case 1
            bitmapData.sensitivity(blockID) = datastruct(currentBlockID).sensitivity;
            bitmapData.adaptthresh(blockID)= NaN;
        case 0
            bitmapData.sensitivity(blockID) = NaN;
            bitmapData.adaptthresh(blockID)=datastruct(currentBlockID).adaptthresh;
    end

    bitmapData.orts(1,:,blockID) = [0 90];
    if ~isempty(datastruct(currentBlockID).gaussianContourLevel) && numel(datastruct(currentBlockID).gaussianContourLevel)<2
        bitmapData.gaussianContourLevel(:,:,blockID)=repmat(datastruct(currentBlockID).gaussianContourLevel,1,2);
    elseif isempty(datastruct(currentBlockID).gaussianContourLevel)
        bitmapData.gaussianContourLevel(:,:,blockID)=repmat(NaN,1,2);
    elseif numel(datastruct(currentBlockID).gaussianContourLevel)==2
        bitmapData.gaussianContourLevel(:,:,blockID)=datastruct(currentBlockID).gaussianContourLevel;
    end
    if numel(datastruct(currentBlockID).gaussianContourLevelMax)<2 && ~isempty(datastruct(currentBlockID).gaussianContourLevelMax)
        bitmapData.gaussianContourLevelMax(:,blockID)=repmat(datastruct(currentBlockID).gaussianContourLevelMax,1,2);
    elseif isempty(datastruct(currentBlockID).gaussianContourLevelMax)
        bitmapData.gaussianContourLevelMax(:,blockID)=repmat(NaN,1,2);
    else
        bitmapData.gaussianContourLevelMax(:,blockID)=datastruct(currentBlockID).gaussianContourLevelMax;
    end
    if isempty(datastruct(currentBlockID).gaussianCond)
        bitmapData.gaussianCond(blockID)=NaN;
    else
        bitmapData.gaussianCond(blockID)=datastruct(currentBlockID).gaussianCond;
    end
    if isempty(datastruct(currentBlockID).nColumnsWanted)
        datastruct(currentBlockID).nColumnsWanted=[NaN NaN];
    end
    bitmapData.nColumnsWanted(:,:,blockID)=datastruct(currentBlockID).nColumnsWanted;
    bitmapData.projectorx=1920;
    bitmapData.projectory=1080;
    bitmapData.pixelsizemm=.0054;
    bitmapData.ProjTTLPulseOn(blockID)=datastruct(currentBlockID).powercycle*50;
    bitmapData.ProjTTLPulseOff(blockID)=50-(datastruct(currentBlockID).powercycle*50);

    
    % copy original bitmaps to run folder
    cloneLoadFlag='load';
    cloneBitmaps(currentBlockStruct, cloneLoadFlag);
end

