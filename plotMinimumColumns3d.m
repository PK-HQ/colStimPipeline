function plotMinimumColumns3d(xColumns,yBeta,zPower,analysisSessID)
% Separate the actual blocks from control for analysis and plotting
controlBlocks=[find(analysisSessID==50) find(analysisSessID==51)]; %identify control blocks
zPowerActual=zPower(setdiff(1:numel(analysisSessID),controlBlocks),:);
zPowerControl=zeros(2,numel(controlBlocks));
xColumnsActual=xColumns(setdiff(1:numel(analysisSessID),controlBlocks),:);
xColumnsControl=zeros(2,numel(controlBlocks));
yBetaActual=yBeta(setdiff(1:numel(analysisSessID),controlBlocks),:);
yBetaControl=yBeta(controlBlocks,:);

%% Fit surface for all points
figure;[hAx,~]=tight_subplot(1,3,[.1 .15],[.1 .1],[.1 .1]);
patternStr={'H','V'};
for pattern=1:2
    axes(hAx(pattern))
    x=rmnan(xColumnsActual(:,pattern));
    y=rmnan(zPowerActual(:,pattern));
    z=rmnan(yBetaActual(:,pattern)*100);
    [surfFit,gof] = fit([x y],z,'poly22','Normalize','off')
    surfPlot=plot(surfFit, [x y],z); hold on
    xlim([0 50]); ylim([0 10]*10^4); zlim([-100 100])
    surfPlot(1).FaceAlpha=.5;
    surfPlot(1).MarkerFaceColor='r';
    surfPlot(1).MarkerEdgeColor=[0 0 0];
    surfPlot(1).MarkerSize=20;
    colormap(plasma)

    ylim([0 22] *10^4)
    title(strcat('Optostim_{', patternStr(pattern), '}'))
    xlabel('Columns');ylabel('Energy');zlabel(strcat('\beta_{', patternStr(pattern), '}'))
%     legend({'optostim'})
    upFontSize(14,.0025)
end
%% Same for controlled power
for pattern=3
    axes(hAx(pattern))
    
    x=rmnan(mean([xColumnsActual(:,2),xColumnsActual(:,1)],2));
    y=rmnan(mean([zPowerActual(:,2),zPowerActual(:,1)],2));
    z=rmnan(yBetaActual(:,2)*100-yBetaActual(:,1)*100);
    
    [surfFit,gof] = fit([x y],z,'poly22','Normalize','off')
    surfPlot=plot(surfFit,[x y],z); hold on
    xlim([0 50]); ylim([0 10]*10^4); zlim([-100 100])
    surfPlot(1).FaceAlpha=.5;
    surfPlot(1).MarkerFaceColor='r';
    surfPlot(1).MarkerEdgeColor=[0 0 0];
    surfPlot(1).MarkerSize=20;
    colormap(plasma)
    ylim([0 22] *10^4)

    title('Columns x biasing effect (fixed power)')
    xlabel('Columns');ylabel('Energy');zlabel('\Delta\beta_{V-H}')
%     legend({'optostim'})
    upFontSize(14,.0025)
end
suplabel('Columns x biasing effect (all powers)','t')



%% Same for controlled power
% Detect blocks with outliers in stimulation energy (0 and 90deg pattern-wise), remove from all xyz
outlierBlocks=isoutlier(zPowerActual,'median'); % 
zPowerActual(outlierBlocks)=NaN;
xColumnsActual(outlierBlocks)=NaN;
yBetaActual(outlierBlocks)=NaN;
zPowerRescale=rescale(zPowerActual,.25,1); %from 0.5-1.0 for plotting

figure;[hAx,~]=tight_subplot(1,3,[.1 .15],[.1 .1],[.1 .1]);
patternStr={'H','V'};
for pattern=1:2
    axes(hAx(pattern))
    
    x=rmnan(xColumnsActual(:,pattern));
    y=rmnan(zPowerActual(:,pattern));
    z=rmnan(yBetaActual(:,pattern)*100);
    
    [surfFit,gof] = fit([x y],z,'poly22','Normalize','off')
    surfPlot=plot(surfFit,[x y],z); hold on
    xlim([0 50]); ylim([0 10]*10^4); zlim([-100 100])
    surfPlot(1).FaceAlpha=.5;
    surfPlot(1).MarkerFaceColor='r';
    surfPlot(1).MarkerEdgeColor=[0 0 0];
    surfPlot(1).MarkerSize=20;
    colormap(plasma)
    ylim([0 22] *10^4)

    title(strcat('Optostim_{', patternStr(pattern), '}'))
    xlabel('Columns');ylabel('Energy');zlabel(strcat('\beta_{', patternStr(pattern), '}'))
%     legend({'optostim'})
    upFontSize(14,.0025)
end

% Same for controlled power
for pattern=3
    axes(hAx(pattern))
    
    x=[mean([xColumnsControl(:,2),xColumnsControl(:,1)],2); rmnan(mean([xColumnsActual(:,2),xColumnsActual(:,1)],2))];
    y=[mean([zPowerControl(:,2),zPowerControl(:,1)],2); rmnan(mean([zPowerActual(:,2),zPowerActual(:,1)],2))];
    z=[yBetaControl(:,2)*100-yBetaControl(:,1)*100; rmnan(yBetaActual(:,2)*100-yBetaActual(:,1)*100)];
    
    
    [mdlfit,gof] = fit([x y],z,'poly22','robust','off')
    surfPlot=plot(mdlfit,[x y],z); hold on
    xlim([0 50]); ylim([0 10] *10^4); zlim([-10 100]); 
    %surfPlot(1).FaceAlpha=.5;
    %surfPlot(1).MarkerFaceColor='r';
    %surfPlot(1).MarkerEdgeColor=[0 0 0];
    %surfPlot(1).MarkerSize=20;
    colormap(plasma)
    %ylim([0 22] *10^4)

    title('Columns x biasing effect (fixed power)')
    xlabel('Columns');ylabel('Power');zlabel('\Delta\beta_{V-H}')
%     legend({'optostim'})
    upFontSize(14,.0025)
end
suplabel('Columns x biasing effect (controlled power)','t')


%% Same for controlled power
for pattern=3
    axes(hAx(pattern))
    
    x=rmnan(mean([xColumnsActual(:,2),xColumnsActual(:,1)],2));
    y=rmnan(mean([zPowerActual(:,2),zPowerActual(:,1)],2));
    z=rmnan(yBetaActual(:,2)*100-yBetaActual(:,1)*100);
    
    [surfFit,gof] = fit([x y],z,'poly22','Normalize','on')
    surfPlot=plot(surfFit,[x y],z); hold on
    xlim([0 50]); ylim([0 10]*10^4); zlim([-100 100])
    surfPlot(1).FaceAlpha=.5;
    surfPlot(1).MarkerFaceColor='r';
    surfPlot(1).MarkerEdgeColor=[0 0 0];
    surfPlot(1).MarkerSize=20;
    colormap(plasma)
    ylim([0 22] *10^4)

    title('Columns x biasing effect (fixed power)')
    xlabel('Columns');ylabel('Energy');zlabel('\Delta\beta_{V-H}')
%     legend({'optostim'})
    upFontSize(14,.0025)
end



%
end
