% binning (bin window=10, centered at 10s)
edgesMin = [-105:10:95];
edgesMax = [-95:10:105];
edgeCenters=mean([edgesMin;edgesMax]);
betaCenters=nan(3,numel(edgeCenters),size(psychometrics,2));
nConds=3;
for nSessions=1:size(psychometrics,2)
    sessionContrast=psychometrics(nSessions).contrastsSorted;
    sessionBeta=psychometrics(nSessions).percentVerticalSorted;

    for ii=1:nConds
        condContrast=sessionContrast(:,ii);
        condBeta=sessionBeta(:,ii);

        % for each bin, get average beta for that bin
        for edge=1:numel(edgesMin)
            edgeMin=edgesMin(edge); edgeMax=edgesMax(edge);
            idx=find(condContrast>=edgeMin & condContrast<edgeMax);
            if ~isempty(idx)
                betaCenters(ii,edge,nSessions)=nanmean(condBeta(idx));
            end
        end
        
    end
end
    
%% Mean and SEM across pages
sessionBinnedContrasts=edgeCenters; % total contrasts, in bins of 10, centered around 10s
sessionBetaMean=nanmean(betaCenters,3); % mean across Sessions
betaCenterSize=nan(nConds,numel(edgesMin));
for ii=1:nConds
    for edge=1:numel(edgesMin)
        betas=squeeze(betaCenters(ii,edge,:));
        betaCenterSize(ii,edge)=numel(betas(~isnan(betas)));
    end
end
betaCenterSize(betaCenterSize==0)=NaN;
SEM = nanstd(betaCenters,[],3)./sqrt(betaCenterSize); % SEM across Sessions

%% Percent vertical vs binned contrasts
col={'b','k','r'};
figure('name','Vertical reports across sessions')
for cond=1:size(sessionBetaMean,1)
    errorbar(sessionBinnedContrasts, sessionBetaMean(cond,:), SEM(cond,:),'Color',col{cond},'LineWidth',3); hold on
end
ylabel('Percent vertical (%)'); xlabel('Binned contrasts (%)')
upFontSize(32,.008); axis square
xline(0,'--','Color',[.4 .4 .4])
yline(50,'--','Color',[.4 .4 .4])
legend({'H-opto','Baseline','V-opto'},'Location','northwest')
title(sprintf('Psychometric data (n_{sessions}=%d)',size(psychometrics,2)),'FontWeight','Normal')

%% Percent horizontal vs binned contrasts
col={'b','k','r'};
figure('name','Vertical reports across sessions')
for cond=1:size(sessionBetaMean,1)
    errorbar(sessionBinnedContrasts, sessionBetaMean(cond,:), SEM(cond,:),'Color',col{cond},'LineWidth',3); hold on
end
ylabel('Percent vertical (%)'); xlabel('Binned contrasts (%)')
upFontSize(32,.008); axis square
xline(0,'--','Color',[.4 .4 .4])
yline(50,'--','Color',[.4 .4 .4])
legend({'H-opto','Baseline','V-opto'},'Location','northwest')
title(sprintf('Psychometric data (n_{sessions}=%d)',size(psychometrics,2)),'FontWeight','Normal')