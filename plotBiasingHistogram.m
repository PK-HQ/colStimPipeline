function plotBiasingHistogram(deltaSorted)
%% Delta histogram
% hardcode        (1-seg significant PMC)                                (2-seg significant PMC)
h_opto=deltaSorted(:,1);
bl_opto=deltaSorted(:,2);
v_opto=deltaSorted(:,3);

%% Run MWU test, calculate p-value that both dist different, 2-tailed
[p_hv,~]=ranksum(h_opto,v_opto,'tail','both')
[p_hbl,~]=ranksum(h_opto,bl_opto,'tail','both')
[p_vbl,~]=ranksum(v_opto,bl_opto,'tail','both')

%% Plot distributions and display medians and MWU p-value
figure
xlim([0 1])
ylim([0 10])

%h
histogram(h_opto,'FaceColor','b','FaceAlpha',.6,'BinWidth',0.05,'LineWidth',2); hold on
%bl
histogram(bl_opto,'FaceColor',[.85 .85 .85],'FaceAlpha',.6,'BinWidth',0.05,'LineWidth',2); hold on
%v
histogram(v_opto,'FaceColor','r','FaceAlpha',.6,'BinWidth',0.05,'LineWidth',2); hold on

xlim([0 1])
ylim([0 8])

%median lines
scatter(median(h_opto),6.8,150,'vk','filled','HandleVisibility','off')
scatter(median(h_opto),6.8,60,'v','filled','MarkerFaceColor','b')
scatter(median(bl_opto),6.8,150,'vk','filled','HandleVisibility','off')
scatter(median(bl_opto),6.8,60,'v','filled','MarkerFaceColor',[.85 .85 .85])
scatter(median(v_opto),6.8,150,'vk','filled','HandleVisibility','off')
scatter(median(v_opto),6.8,60,'vk','filled','MarkerFaceColor','r')


xlabel('\Delta')
ylabel('Count')
addSkippedTicks(0,1,.1,'x')

addSkippedTicks(0,8,1,'y')
hLegend=legend({'H-opto','Baseline','V-opto',...
    },'location','eastoutside');
hTitle=title('\Delta distributions for stimulation conditions','FontWeight','Normal');

upFontSize(16,0.01);
hLegend.FontSize=14;
text((median(h_opto)+median(v_opto))/2, 7.4, '***','FontSize',20)

text((median(h_opto)+median(bl_opto))/2, 7.8, 'n.s.','FontSize',20)

text((median(v_opto)+median(bl_opto))/2, 7.4, '***','FontSize',20)
end