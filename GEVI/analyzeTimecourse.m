function analyzeTimecourse(dataStruct)
    %set(0,'DefaultFigureWindowStyle','docked')
    switch dataStruct.series
        case {'Contrast'}
            legStr={'100% (0\circ)', '50%', '25%', '12.5%', '6.25%', '3.125%',...
                    '100% (90\circ)', '50%', '25%', '12.5%', '6.25%', '3.125%'};
        case {'Temporal'}
            legStr={'2Hz (0\circ)', '4Hz', '8.3Hz', '10Hz', '16.7Hz', '20Hz', '33.3Hz',...
                    '2Hz (90\circ)', '4Hz', '8.3Hz', '10Hz', '16.7Hz', '20Hz', '33.3Hz'};
    end   
    % load
    filenameStruct=generateFilenames(dataStruct);
    load(filenameStruct.vdaq)
    load(filenameStruct.ROITC)
    load(filenameStruct.TS)

    
    % def ROI
    xROI=Crop(1):Crop(1)+Crop(3);
    yROI=Crop(2):Crop(2)+Crop(4);

    %def blank and resp, and for ROI
    blankConds=find(TS.Header.Conditions.TypeCond==0);
    stimConds=find(TS.Header.Conditions.TypeCond==2);

    nStimConds=numel(stimConds);
    
    blankResp=mean(RespCond(:,:,blankConds),3);
    condResp=RespCond(:,:,stimConds)-blankResp;
    condRespTest=condResp;
    blankRespROI=blankResp(yROI,xROI);
    condRespROI=condResp(yROI,xROI,:);

    %% Plot ROI
    figure('name','Raw cond. responses')
    [hAx,~]=tight_subplot(2,nStimConds/2,[.015 .015]);%delete(hAx([4 12]));
    for i=1:size(condResp,3)
        axes(hAx(i))
        imgsc(condResp(:,:,i),legStr{i}); colorbar; hold on
        rectangle('Position',[Crop(1), Crop(2), Crop(3), Crop(4)],...
             'LineWidth',2,'LineStyle','-','EdgeColor','r')
    end
    suplabel([dataStruct.series ' series: Raw ROI image'],'t')

    %% Plot ROI TC
    condTC=DataCond(:,:,:,stimConds)-mean(DataCond(:,:,:,blankConds),4);
    figure('name','Raw cond. timecourses')
    
    % def n stimconds
    nColors=nStimConds+4;
    nStimConds=nStimConds;
    
    cmap=redblue(nStimConds+6);
    colormap(1:nStimConds/2,:)=cmap(1:nStimConds/2,:);
    colormap(1+nStimConds/2:nStimConds,:)=flipud(cmap(2+length(cmap)/2:1+length(cmap)/2+nStimConds/2,:));

    for i=1:size(condTC,4)
        plot(1:140,squeeze(mean(condTC(yROI,xROI,:,i),[1 2]))*100,'LineWidth',2,'LineStyle','-','Color',colormap(i,:)); hold on
    end
    ylabel('\Delta F/F (%)');addSkippedTicks(0,140,10,'x')
    xlabel('Time (s)');addSkippedTicks(0,140,10,'x')
    suplabel([dataStruct.series ' series: Raw ROI timecourse'],'t')
    upFontSize(20,0.005)
    hLeg=legend(legStr,'NumColumns',2,'FontSize',12,'position','northwest');
    title(hLeg,dataStruct.series)
    
    %% Plot ROI TC averaged
    
    figure('name','Avg cond. timecourses')
    
    % average first and second pulse
    stimEpoch=1000;
    stimWavelength=TS.Header.Conditions.DurationOn+TS.Header.Conditions.DurationOff;
    stimHz=stimEpoch./stimWavelength;
    daqHz=round(VDaqSettings.datalog.framerate);
    exposureTime=round(VDaqSettings.datalog.expotime);
    nFrames=VDaqSettings.cam.vid.FramesPerTrigger;
    
    nmlFrame=13; %13 if 100Hz
    startFrame=nmlFrame+1; %14
    endFrame=nFrames; %140 or 120
    nCycles=stimHz; % >=2Hz
    stimWavelengthFrame=stimWavelength./exposureTime; % 50 if 2Hz, 25 if 4Hz ...
    
    
    % change frames to non hardcoded
    condTCAvg=nan(size(condTC,1),size(condTC,2),100,nStimConds);
    for stimCond=1:nStimConds
        condTCAvg(:,:,1:1+stimWavelengthFrame(stimCond)-1,stimCond)=...
            (condTC(:,:,startFrame:startFrame+stimWavelengthFrame(stimCond)-1,stimCond)+condTC(:,:,startFrame+stimWavelengthFrame(stimCond):startFrame+2*stimWavelengthFrame(stimCond)-1,stimCond))./2;
    end

    for i=6:12%nStimConds
        x=1:100;
        y=squeeze(nanmean(condTCAvg(yROI,xROI,:,i),[1 2]))*100;
        %idx = ~isnan(y);
        %y=y(idx);
        plot(x,y,'LineWidth',2,'LineStyle','-','Color',colormap(i,:)); hold on
    end
    ylabel('Average \Delta F/F (%)');addSkippedTicks(0,140,10,'x')
    xlabel('Time (s)');addSkippedTicks(0,140,10,'x')
    suplabel([dataStruct.series ' series: Averaged ROI timecourse'],'t')
    upFontSize(20,0.005)
    hLeg=legend(legStr,'NumColumns',2,'FontSize',14);
    title(hLeg,dataStruct.series)

    
end