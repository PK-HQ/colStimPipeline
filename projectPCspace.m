function projectedDataCont = projectPCspace(corticalResponseTrial,trialConditions,blankAveraged,nOptoConds,nBlank,frames,PCACoef,nPCAComponents)
topPC=PCACoef(:,1:nPCAComponents);
for optoCondID=1:nOptoConds
    [condTrialIdx,~]=find(trialConditions==optoCondID+nBlank); %to get opto cond ID, use 1:nOptoConds + number of blank
    % get condition response averaged across all trials, blank subtracted (pixX, pixY, nFrames, nTrials)
    condResp=mean(corticalResponseTrial(:,:,:,condTrialIdx)-blankAveraged,4);
    counter=1;
    for frame=frames
        % reshape 64 px x 64 px x 1 frame x 317 trials into 4096 x 317 trials
        reshapedData=reshape(condResp(:,:,frame),...
        [size(condResp(:,:,frame),1)*size(condResp(:,:,frame),2),size(condResp(:,:,frame),3)])';
        for nPC=1:nPCAComponents
          % project data into 3-dimensional PC space (top-3 components)
          projectedData = reshapedData*topPC(:,nPC);
          % store as frames x opto cond x PC-coordinates 1-3
          projectedDataCont(counter,optoCondID,nPC)=projectedData;
        end
        counter=counter+1;
    end
end
