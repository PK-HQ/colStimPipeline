function plotMinimumColumns(xColumns,yBeta,zPower,analysisSessID,metricStr)
% Separate the actual blocks from control for analysis and plotting
controlBlocks=[find(analysisSessID==50) find(analysisSessID==51)]; %identify control blocks
zPowerActual=zPower(setdiff(1:numel(analysisSessID),controlBlocks),:);
zPowerControl=zeros(2,numel(controlBlocks));
xColumnsActual=xColumns(setdiff(1:numel(analysisSessID),controlBlocks),:);
xColumnsControl=zeros(2,numel(controlBlocks));
yBetaActual=yBeta(setdiff(1:numel(analysisSessID),controlBlocks),:);
yBetaControl=yBeta(controlBlocks,:);

% Detect blocks with outliers in stimulation energy (0 and 90deg pattern-wise), remove from all xyz
[outlierBlocks,minMAD,maxMAD,centerMAD]=isoutlier(zPowerActual); % 
zPowerActual(outlierBlocks)=NaN;
xColumnsActual(outlierBlocks)=NaN;
yBetaActual(outlierBlocks)=NaN;

%rm #2
zPowerActual(2)=NaN;
xColumnsActual(2)=NaN;
yBetaActual(2)=NaN;

zPowerRescale=rescale(zPowerActual,.15,1); %from 0.5-1.0 for plotting



% Titles, text and legend
switch metricStr
    case {'beta'}
        
        % Plot for H and V
        figure('name','Minimum columns (\beta)');hold on        
        yline(50,'--','HandleVisibility','off'); hold on % reference line for no biasing
        colors={'blue','red'};
        for pattern=1:size(zPowerActual,2)
            color=colors{pattern};
            legendFlag=1;
            for block=1:length(xColumnsActual(:,pattern))
                if ~isnan(zPowerRescale(block,pattern))
                    if legendFlag==1
                        hScatter=scatter(xColumnsActual(block,pattern),yBetaActual(block,pattern)*100,250,color,'filled','MarkerFaceAlpha',zPowerRescale(block,pattern),'MarkerEdgeColor','k','LineWidth',2); hold on;
                        legendFlag=0;
                    else
                        hScatter=scatter(xColumnsActual(block,pattern),yBetaActual(block,pattern)*100,250,color,'filled','MarkerFaceAlpha',zPowerRescale(block,pattern),'MarkerEdgeColor','k','LineWidth',2,'HandleVisibility','off'); hold on;
                    end
                end
            end
            %scatter(xColumnsControl(:,pattern),yBetaControl(:,pattern)*100,150,color,'filled','x','MarkerFaceAlpha',1,'MarkerEdgeColor',color,'LineWidth',2,'HandleVisibility','off'); hold on;
        end
        hold off
        
        % Titles, text, legend
        title('Columns x biasing effect (\beta)')
        xlabel('No. of columns per bitmap','FontName','Arial')
        legend({'H-optostim','V-optostim'})
        upFontSize(18,.0025)

        ylabel('\beta','FontName','Arial','FontSize',24,'FontWeight','bold')
        ylim([0 1]*100); yticks([0:.1:1]*100); 
        xlim([0 40])
    case {'deltabeta'}
        
        % compute mean columns, power and \delta\beta
        zPowerMean=mean(zPowerActual,2);
        xColumnsMean=mean(xColumnsActual,2);
        yDeltaBeta=yBetaActual(:,2)-yBetaActual(:,1);
        
        %yDeltaBetaControl=yBetaControl(:,2)-yBetaControl(:,1);
        %xControl=xColumnsControl(:,2)-xColumnsControl(:,1);
        %zControl=zPowerControl(:,2)-zPowerControl(:,1);

        zPowerRescaleMean=rescale(zPowerMean,.25,1); %from 0.5-1.0 for plotting
        
        % Plot for H and V
        figure('name','Minimum columns (\Delta\beta)');hold on
        yline(0,'--','HandleVisibility','off'); hold on % reference line for no biasing
        %colors={rgb(138,43,226)};
        for pattern=1:size(zPowerMean,2)
            color=[255, 0, 153]/255;
            legendFlag=1;
            for block=1:length(xColumnsMean(:,pattern))
                if ~isnan(zPowerRescaleMean(block,pattern))
                    alphaVal=1;%zPowerRescaleMean(block,pattern);
                    if legendFlag==1
                        hScatter=scatter(xColumnsMean(block,pattern),yDeltaBeta(block,pattern),250,color,'filled','^','MarkerFaceAlpha',alphaVal,'MarkerEdgeColor','k','LineWidth',2); hold on;
                        legendFlag=0;
                    else
                        hScatter=scatter(xColumnsMean(block,pattern),yDeltaBeta(block,pattern),250,color,'filled','^','MarkerFaceAlpha',alphaVal,'MarkerEdgeColor','k','LineWidth',2,'HandleVisibility','off'); hold on;
                    end
                end
            end
            %scatter(xColumnsControl(:,pattern),yDeltaBetaControl(:,pattern)*100,150,color,'filled','x','MarkerFaceAlpha',1,'MarkerEdgeColor',color,'LineWidth',2,'HandleVisibility','off'); hold on;
        end
        hold off
        
        title('Columns x biasing effect (\Delta\beta)','FontWeight','Normal')
        xlabel('Average no. of columns across bitmaps','FontName','Arial')
        legend({'\Delta\beta'},'FontWeight','bold')
        upFontSize(18,.0025)

        ylabel('\Delta\beta','FontSize',24,'FontWeight','bold')
        ylim([-.05 .6]); yticks([-.1:.1:.6]); 
        xlim([0 40]) 
    case {'deltaphi'}
        
        % compute mean columns, power and \delta\beta
        zPowerMean=mean(zPowerActual,2);
        xColumnsMean=mean(xColumnsActual,2);
        yDeltaBeta=yBetaActual(:,2)-yBetaActual(:,1);
        
        %yDeltaBetaControl=yBetaControl(:,2)-yBetaControl(:,1);
        %xControl=xColumnsControl(:,2)-xColumnsControl(:,1);
        %zControl=zPowerControl(:,2)-zPowerControl(:,1);

        zPowerRescaleMean=rescale(zPowerMean,.25,1); %from 0.5-1.0 for plotting
        
        % Plot for H and V
        figure('name','Minimum columns (\Delta\phi)');hold on
        yline(0,'--','HandleVisibility','off'); hold on % reference line for no biasing
        %colors={rgb(138,43,226)};
        for pattern=1:size(zPowerMean,2)
            color=[255, 192, 0]/255;
            legendFlag=1;
            for block=1:length(xColumnsMean(:,pattern))
                if ~isnan(zPowerRescaleMean(block,pattern))
                    alphaVal=1;%zPowerRescaleMean(block,pattern);
                    if legendFlag==1
                        hScatter=scatter(xColumnsMean(block,pattern),yDeltaBeta(block,pattern),250,color,'filled','^','MarkerFaceAlpha',alphaVal,'MarkerEdgeColor','k','LineWidth',2); hold on;
                        legendFlag=0;
                    else
                        hScatter=scatter(xColumnsMean(block,pattern),yDeltaBeta(block,pattern),250,color,'filled','^','MarkerFaceAlpha',alphaVal,'MarkerEdgeColor','k','LineWidth',2,'HandleVisibility','off'); hold on;
                    end
                end
            end
            %scatter(xColumnsControl(:,pattern),yDeltaBetaControl(:,pattern)*100,150,color,'filled','x','MarkerFaceAlpha',1,'MarkerEdgeColor',color,'LineWidth',2,'HandleVisibility','off'); hold on;
        end
        hold off
        
        title('Columns x biasing effect (\Delta\beta)','FontWeight','Normal')
        xlabel('Average no. of columns across bitmaps','FontName','Arial')
        legend({'\Delta\phi'},'FontWeight','bold')
        upFontSize(18,.0025)

        ylabel('\Delta\phi','FontSize',24,'FontWeight','bold')
        ylim([-1 6]); yticks([-1:1:6]); 
        xlim([0 40]) 
        
    case {'deltabetaphi'}
        
        % compute mean columns, power and \delta\beta
        zPowerMean=mean(zPowerActual,2);
        xDeltaPhi=xColumnsActual(:,2)-xColumnsActual(:,1);
        yDeltaBeta=yBetaActual(:,2)-yBetaActual(:,1);
        
        %yDeltaBetaControl=yBetaControl(:,2)-yBetaControl(:,1);
        %xControl=xColumnsControl(:,2)-xColumnsControl(:,1);
        %zControl=zPowerControl(:,2)-zPowerControl(:,1);

        zPowerRescaleMean=rescale(zPowerMean,.25,1); %from 0.5-1.0 for plotting
        
        % Plot for H and V
        figure('name','Minimum columns (\Delta\phi)');hold on
        yline(0,'--','HandleVisibility','off'); hold on % reference line for no biasing
        %colors={rgb(138,43,226)};
        for pattern=1:size(zPowerMean,2)
            color='green';
            legendFlag=1;
            for block=1:length(xDeltaPhi(:,pattern))
                if ~isnan(zPowerRescaleMean(block,pattern))
                    alphaVal=1;%zPowerRescaleMean(block,pattern);
                    if legendFlag==1
                        hScatter=scatter(xDeltaPhi(block,pattern),yDeltaBeta(block,pattern),250,color,'filled','^','MarkerFaceAlpha',alphaVal,'MarkerEdgeColor','k','LineWidth',2); hold on;
                        legendFlag=0;
                    else
                        hScatter=scatter(xDeltaPhi(block,pattern),yDeltaBeta(block,pattern),250,color,'filled','^','MarkerFaceAlpha',alphaVal,'MarkerEdgeColor','k','LineWidth',2,'HandleVisibility','off'); hold on;
                    end
                end
            end
            %scatter(xColumnsControl(:,pattern),yDeltaBetaControl(:,pattern)*100,150,color,'filled','x','MarkerFaceAlpha',1,'MarkerEdgeColor',color,'LineWidth',2,'HandleVisibility','off'); hold on;
        end
        hold off
        
        title('\Delta\phi x \Delta\beta','FontWeight','Normal')
        xlabel('\Delta\phi','FontName','Arial')
       
        upFontSize(18,.0025)

        ylabel('\Delta\beta','FontSize',24,'FontWeight','bold')
        ylim([0 .8]); addSkippedTicks(0,0.8,.1,'y') 
        xlim([0 6]); addSkippedTicks(0,6,1,'x') 

end

end


