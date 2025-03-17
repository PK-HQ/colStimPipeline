function mdlStruct=analyzePsychometricModels(monkeyName, chamberWanted, modelTypes, mainPath,...
    behavioralData, bitmapData, datastruct, analysisBlockID, clusterIdx, plotFlag, saveFlag)
%% Plot by cluster, and by model type
nClusters=sort(clusterIdx,'descend');%unique(clusterIdx);
% Initialize mdlStruct and other structures
mdlStruct = struct();
AICCdeltax = struct();
AICCbeta = struct();
thresh=[];beta=[];exp=[];c50=[];
for cluster=3%1:nClusters%1:nClusters
    disp(['Cluster ' num2str(cluster)])
    
    % Get columns
    columnsDesired=20; columnSpread=4;
    
    % Compute columns mean
    columns = mean(bitmapData.nColumns', 2)';

    % Compute bitmap energy mean
    bitmapEnergy = mean(squeeze(bitmapData.energy), 1);

    % Find cluster block indices matching criteria
    clusterBlocksIdx = find(clusterIdx == cluster); %find(clusterIdx == cluster & ...
        %columns >= (columnsDesired - columnSpread) & ...
        %columns <= (columnsDesired + columnSpread));
    nClusterBlocks=numel(clusterBlocksIdx);

    if nClusterBlocks==0
        continue
    else
        columns=columns(clusterBlocksIdx);
        bitmapEnergy=bitmapEnergy(clusterBlocksIdx);
    end
    bitmapDataAvg.energy=[];
    
    for modelID = 4%8
        modelTypeStr = modelTypes{modelID};

        conds = 1:3;
        xBlocks = squeeze(behavioralData.gaborContrasts(conds, :, clusterBlocksIdx));
        yBlocks = squeeze(behavioralData.percentageCorrect(conds, :, clusterBlocksIdx));
        
        % Fit psychometric data
        switch strcmp(modelTypeStr, 'bill')
            case 1
                [mdl, mdlAvg] = fitBayesianModelMLE(xBlocks, yBlocks, modelTypeStr);
            case 0
                [mdl,mdlAvg] = fitNakaRushtonMLE3(xBlocks, yBlocks, modelTypeStr);
        end
        xFit = sort([-nlinspace(0, 100, 100, 'nonlinear'), nlinspace(0, 100, 100, 'nonlinear')]);
        savefilenameBlock = ['psychometrics/' chamberWanted '-chamber/' modelTypeStr '/' 'psychfit-' chamberWanted '-' modelTypeStr '-C' num2str(cluster)];
        plotFlag=1;
        switch plotFlag
            case 1
                if modelID>0
                    % Plot per block
                    plotAverageFlag=0;
                    mdl=plotNakaRushtonFit4(behavioralData, bitmapData, datastruct, analysisBlockID,...
                        mdl, mdl.fittedParams(:, :, 1), xFit, monkeyName, clusterBlocksIdx, plotAverageFlag,...
                        saveFlag, cluster, modelTypeStr, savefilenameBlock);
                    
                    % Plot block average
                    %{
                    savefilenameAverage = ['psychometrics/' chamberWanted '-chamber/' modelTypeStr '/' 'psychfit-' chamberWanted '-' modelTypeStr '-C' num2str(cluster) '-' num2str(nClusterBlocks+1)];
                    plotAverageFlag=1;
                    plotNakaRushtonFit4(behavioralData, bitmapData, datastruct, analysisBlockID,...
                        mdlAvg, (mdlAvg.fittedParams(:, :, 1)), xFit, monkeyName, clusterBlocksIdx, plotAverageFlag, ...
                        saveFlag, cluster, modelTypeStr, savefilenameAverage);
                    
                    % Save to pdf
                    append_pdfs([mainPath monkeyName '/Meta/' savefilenameBlock '.pdf'],...
                        [mainPath monkeyName '/Meta/' savefilenameAverage '.pdf'], [mainPath monkeyName '/Meta/' savefilenameBlock '.pdf'])        
                    delete([mainPath monkeyName '/Meta/' savefilenameAverage '.pdf'])
                    %CF;

                    %}
                    % plot fitting params
                    thresh=[thresh;mdl.thresholdContrast];
                    beta=[beta;mdl.fittedParams(:,1)];
                    exp=[exp;mdl.fittedParams(:,2)];
                    c50=[c50;mdl.fittedParams(:,3)];
                
                elseif modelID == 0
                    %{
                    % Process 'bill' model
                    plotAverageFlag = 0;
                    mdl = plotNakaRushtonFit3_Bill(behavioralData, bitmapData, datastruct, analysisBlockID,...
                        mdl, mdl.fittedParams(:, :, 1), xFit, monkeyName, clusterBlocks, plotAverageFlag,...
                        saveFlag, cluster, modelTypeStr, savefilenameBlock);

                    % Plot block average for 'bill'
                    savefilenameAverage = ['psychometrics/' chamberWanted '-chamber/' modelTypeStr '/' 'psychfit-' chamberWanted '-' modelTypeStr '-C' num2str(cluster) '-' num2str(nClusterBlocks + 1)];
                    plotAverageFlag = 1;
                    plotNakaRushtonFit3_Bill(behavioralData, bitmapData, datastruct, analysisBlockID,...
                        mdlAvg, mdlAvg.fittedParams(:, :, 1), xFit, monkeyName, clusterBlocks, plotAverageFlag,...
                        saveFlag, cluster, modelTypeStr, savefilenameAverage);

                    % Save to PDF
                    append_pdfs([mainPath monkeyName '/Meta/' savefilenameBlock '.pdf'], ...
                        [mainPath monkeyName '/Meta/' savefilenameAverage '.pdf'], [mainPath monkeyName '/Meta/' savefilenameBlock '.pdf']);
                    delete([mainPath monkeyName '/Meta/' savefilenameAverage '.pdf']);
                    %CF;

                    % Collect fitting parameters for summary plots
                    thresh = [thresh; mdl.thresholdContrast];
                    %}
                end

                    
        end
    end
    mdlStruct.([chamberWanted, modelTypeStr, 'C' num2str(cluster)]) = mdl;
    mdlStruct.([chamberWanted, modelTypeStr, 'C' num2str(cluster) 'Avg']) = mdlAvg;
    %% Compare models of given cluster
    % Filter and compare models with non-parametric t-test (per chamber)
    savefilename = ['psychometrics/' chamberWanted '-chamber/' 'paramdist-' chamberWanted '-' modelTypes{2} modelTypes{3}  '-C' num2str(cluster)];
    %bitmapEnergy.(chamberWanted)=squeeze(bitmapData.energy(:,:,clusterBlocksIdx))';
    %[powerFiltered, powerFilteredMask] = powerFiltering(squeeze(bitmapData.energy)');

    %AICCdeltax.(chamberWanted) = mdlStruct.([chamberWanted, 'deltax' , 'C' num2str(cluster)]).fittedParams(:, end);
    %AICCbeta.(chamberWanted) = mdlStruct.([chamberWanted, 'beta', 'C' num2str(cluster)]).fittedParams(:, end);
    %plotParamDistributions(AICCbeta.(chamberWanted), AICCdeltax.(chamberWanted), 'AICc', '/beta vs /Deltax', {'/beta','/Deltax'}, monkeyName, savefilename, saveFlag);    
    
    % Save fitting parameters
    %save([mainPath '/Chip/Meta/psychometrics/fittingParams' chamberWanted '-' modelTypes{2} modelTypes{3}  '-C' num2str(cluster) '.mat'], 'mdlStruct', 'bitmapEnergy', 'AICCdeltax', 'AICCbeta');
    
    % Combined chamber analysis
    %savefilename = [mainPath 'psychometrics/paramdist-' chamberWanted '-' modelTypes{2} modelTypes{3}  '-C' num2str(cluster)];
    %plotParamDistributions([AICCbeta.('L')], [AICCdeltax.('L')], 'AICc', '/beta vs /Deltax', {'/beta','/Deltax'}, monkeyName, savefilename, saveFlag);
    %ranksum(mdlStruct.Rbeta.fittedParams(:, 1), mdlStruct.Lbeta.fittedParams(:, 1));
    %ranksum(mdlStruct.Rdeltax.fittedParams(:, 4), mdlStruct.Ldeltax.fittedParams(:, 4));
    %CF
end
save([mainPath '/Chip/Meta/psychometrics/' chamberWanted '-chamber/' modelTypeStr '/mdlStruct' chamberWanted num2str(columnsDesired) '.mat'], 'mdlStruct');

%% Plot param dist

%{
figure
nRows=1;nCols=4;
gap=.05;marginV=.01;marginH=.075;
[hAx,~]=tight_subplot(nRows,nCols,[gap gap], [marginV+.27 marginV+.13], [marginH marginH]);
axes(hAx(1));
data = unique(thresh);
% Define x-axis bins from 0 to 100 with bin width of 1
% Get counts for each bin
edges = -.5:2.5:50.5;
[counts, binEdges] = histcounts(data, edges);
% Compute bin centers for plotting
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
% Plot the counts
bar(binCenters, counts, 'FaceColor', 'r');
xlabel('Contrast (%)')
ylabel('Counts');
title(['threshold (' num2str(median(data),2) '±' num2str(mad(data),1) ')']);
xlim([0 50]); xticks([0:10:50])
upFontSize(32,.01)
axes(hAx(2));
data = unique(beta);
% Define x-axis bins from 0 to 100 with bin width of 1
% Get counts for each bin
edges = -.5:5:100.5;
[counts, binEdges] = histcounts(data, edges);
% Compute bin centers for plotting
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
% Plot the counts
bar(binCenters, counts, 'FaceColor', 'k');
xlabel('\beta');
title(['\beta (' num2str(median(data),2) '±' num2str(mad(data),1) ')']);
xlim([0 100]); xticks([0:25:100])
upFontSize(32,.01)
axes(hAx(3));
data = unique(exp);
% Define x-axis bins from 0 to 100 with bin width of 1
% Get counts for each bin
edges = -.5:1:6.5;
[counts, binEdges] = histcounts(data, edges);
% Compute bin centers for plotting
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
% Plot the counts
bar(binCenters, counts, 'FaceColor', 'k');
xlabel('n');
title(['exponent (' num2str(median(data),2) '±' num2str(mad(data),1) ')']);
xlim([0 7]); xticks([0:1:7])
upFontSize(32,.01)
axes(hAx(4));
data = unique(c50);
% Define x-axis bins from 0 to 100 with bin width of 1
% Get counts for each bin
edges = -.5:2.5:50.5;
[counts, binEdges] = histcounts(data, edges);
% Compute bin centers for plotting
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
% Plot the counts
bar(binCenters, counts, 'FaceColor', 'k');
xlabel('C50');
title(['C50 (' num2str(median(data),2) '±' num2str(mad(data),1) ')']);
xlim([0 50]); xticks([0:10:50])
upFontSize(32,.01)
suplabel(['L-Chip'],'t',[.1 .1 .85 .85]);
upFontSize(32,.01)
%}
%% Test model fits for significance, per cluster
%{
for cluster=nClusters
    mdl1AIC=mdlStruct.([chamberWanted, modelTypes{2}, 'C' num2str(cluster)]).fittedParams(:,end);
    mdl2AIC=mdlStruct.([chamberWanted, modelTypes{3}, 'C' num2str(cluster)]).fittedParams(:,end);
    k1=sum(~isnan(mdlStruct.([chamberWanted, modelTypes{2}, 'C' num2str(cluster)]).fittedParams(1,:)));
    k2=sum(~isnan(mdlStruct.([chamberWanted, modelTypes{3}, 'C' num2str(cluster)]).fittedParams(1,:)));
    mdlStruct.([chamberWanted 'C' num2str(cluster)]).paramSigTest(cluster)=ranksum(mdl1AIC, mdl2AIC);%ranksum(mdl1Cluster,mdl2Cluster,'tail','right')
    %plotParamDistributions(mdl1AIC, mdl2AIC, 'AICc', 'Beta vs Deltax', {'Beta','Deltax'}, monkeyName, savefilename, saveFlag);
    fprintf('Medians = %.0f & %.0f (p=%.2f)\n', median(mdl1AIC), median(mdl2AIC), ...
        mdlStruct.([chamberWanted 'C' num2str(cluster)]).paramSigTest(cluster))
end
%}
end