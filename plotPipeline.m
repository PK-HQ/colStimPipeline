%% Change these for experiment runs
monkeyName='Pepper';%Pepper or blank
publicationPlot='power-biasing';
% Saving and plotting flags
saveFlag=0;
saveFlagBMP=0;
plotFlag=1;
skipImaging=1;

%% Load dataStruct for the desired chamber
[mainPath, datastruct]=setupEnv(['users/PK/colStimPipeline/exptListBiasingFull' monkeyName '.m']);
chambers={'R', 'L'};
for chamberID=1
    nColumnsWanted=[]; chamberWanted=chambers{chamberID};
    analysisBlockID = organizeBlocks(datastruct, chamberWanted, nColumnsWanted);

    switch publicationPlot
        case {'power-biasing'}
            load([mainPath '/' monkeyName '/Meta/psychometrics/psychfit' chamberWanted '-' modelTypes{2} '.mat'])
            %% 20 column power x biasing
            figure
            conY = mdlStruct.([chamberWanted, 'weibullfreeAll' , 'C1']).yConOpto;
            inconY = mdlStruct.([chamberWanted, 'weibullfreeAll' , 'C1']).yInconOpto;
            deltaY = nanmean(conY-inconY,2);
            xData=mean(squeeze(bitmapData.energy(:,:,:))',2);
            columns=mean(bitmapData.nColumns(:,:))';
            columnsDesired=20; columnSpread=4;
            blocksActual = find(columns >= columnsDesired-columnSpread &...
                columns <= columnsDesired+columnSpread &...
                mean(squeeze(bitmapData.orts)==[0;90])');
            
            blocksControl = find(columns >= columnsDesired-columnSpread &...
                columns <= columnsDesired+columnSpread &...
                mean(squeeze(bitmapData.orts)==[45;135])');
            
            yline(0,'LineStyle','--','color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off'); hold on
            % Plot actual data
            actualX=xData(blocksActual);
            actualY=deltaY(blocksActual);
            mdlScatter=fitSaturatingCurve(actualX, actualY,  'k', 1); hold on
            scatter(actualX,actualY,200,'Marker', 'square', ...
                'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth',2); hold on                
            % Plot control data
            controlX=xData(blocksControl);
            controlY=deltaY(blocksControl);
            if ~isempty(blocksControl)
                scatter(controlX+randi([-10 10],1,numel(controlX))'/100,controlY,200,'Marker', 'square', ...
                    'MarkerFaceColor', [1 1 1]*.7, 'MarkerEdgeColor', 'k', 'LineWidth',2); hold on
            end            
    
            title({'Effect of power on Δcon-incon', sprintf('(%.0f ± %.0f columns)',columnsDesired,round(std(columns)))})
            axis square
            xlim([0 5])
            ylim([-10 50])
            addSkippedTicks(0,5,.5,'x')
            addSkippedTicks(-10,50,5,'y')
            %[p,h,stats] = ranksum(actualY(actualX>=3.5), [controlY(controlX>=3.5)' 2 -1 -3], 'tail','both')
            xlabel('Total energy (mW)')
            ylabel('Δ_{con-incon opto}')
            upFontSize(20,.02)
            %export_fig(['Y:\users\PK\colStimPipeline\figures\powerxdeltaY-' monkeyName '-' chamberWanted],'-svg','-png','-nocrop','-r600');
    
        case {'coreg-biasing'}
            load([mainPath '/' monkeyName '/Meta/psychometrics/psychfit' chamberWanted '-' modelTypes{2} '.mat'])
            %% 20 column power x biasing
            figure
            conY = mdlStruct.([chamberWanted, 'weibullfreeAll' , 'C1']).yConOpto;
            inconY = mdlStruct.([chamberWanted, 'weibullfreeAll' , 'C1']).yInconOpto;
            deltaY = nanmean(conY-inconY,2); 
            xData=bitmapData.similarity;
            columns=mean(bitmapData.nColumns(:,:))';
            columnsDesired=20; columnSpread=4;
            blocksActual = find(columns >= columnsDesired-columnSpread &...
                columns <= columnsDesired+columnSpread &...
                mean(squeeze(bitmapData.orts)==[0;90])');
            
            blocksControl = find(columns >= columnsDesired-columnSpread &...
                columns <= columnsDesired+columnSpread &...
                mean(squeeze(bitmapData.orts)==[45;135])');
            
            yline(0,'LineStyle','--','color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off'); hold on
            % Plot actual data
            actualX=xData(blocksActual);
            actualY=deltaY(blocksActual);
            mdlScatter=fitSaturatingCurve(actualX, actualY,  'k', 1); hold on
            scatter(actualX,actualY,200,'Marker', 'square', ...
                'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth',2); hold on                
            % Plot control data
            controlX=xData(blocksControl);
            controlY=deltaY(blocksControl);
            if ~isempty(blocksControl)
                scatter(controlX+randi([-10 10],1,numel(controlX))/100,controlY,200,'Marker', 'square', ...
                    'MarkerFaceColor', [1 1 1]*.7, 'MarkerEdgeColor', 'k', 'LineWidth',2); hold on
            end            
    
            title({'Effect of coregistration quality on Δcon-incon', sprintf('(%.0f ± %.0f columns)',columnsDesired,round(std(columns)))})
            axis square
            xlim([0 1])
            ylim([-10 50])
            addSkippedTicks(0,1,.25,'x')
            addSkippedTicks(-10,50,5,'y')
            %[p,h,stats] = ranksum(actualY(actualX>=3.5), [controlY(controlX>=3.5)' 2 -1 -3], 'tail','both')
            xlabel('SSIM (a.u.)')
            ylabel('Δ_{con-incon opto}')
            upFontSize(20,.02)
    end
end