function plotBiasingHistogram(mainPath,psychometric,saveFlag)
%% Delta histogram
h_opto=psychometric.betaSortedCont(:,1);
bl_opto=psychometric.betaSortedCont(:,2);
v_opto=psychometric.betaSortedCont(:,3);

%% Run MWU test, calculate p-value that both dist different, 2-tailed
[p_hv,~]=ranksum(h_opto,v_opto,'tail','both')
[p_hbl,~]=ranksum(h_opto,bl_opto,'tail','both')
[p_vbl,~]=ranksum(v_opto,bl_opto,'tail','both')

%% Plot distributions and display medians and MWU p-value
figure
xlim([0 1])
ylim([0 10])

%h-opto histogram
histogram(h_opto,'FaceColor','b','FaceAlpha',.6,'BinWidth',0.05,'LineWidth',2); hold on
scatter(median(h_opto),5.7,150,'vk','filled','HandleVisibility','off') %median mkr
scatter(median(h_opto),5.7,60,'v','filled','MarkerFaceColor','b','HandleVisibility','off')


%v-opto histogram
histogram(v_opto,'FaceColor','r','FaceAlpha',.6,'BinWidth',0.05,'LineWidth',2); hold on
scatter(median(v_opto),5.7,150,'vk','filled','HandleVisibility','off')
scatter(median(v_opto),5.7,60,'vk','filled','MarkerFaceColor','r','HandleVisibility','off')
xlim([0 1])
ylim([0 6])

%bl-opto histogram
%{
histogram(bl_opto,'FaceColor',[.85 .85 .85],'FaceAlpha',.6,'BinWidth',0.05,'LineWidth',2); hold on
scatter(median(bl_opto),6.8,150,'vk','filled','HandleVisibility','off')
scatter(median(bl_opto),6.8,60,'v','filled','MarkerFaceColor',[.85 .85 .85])
%}

xlabel('\beta')
ylabel('Count')
addSkippedTicks(0,1,.1,'x')

addSkippedTicks(0,8,1,'y')
hLegend=legend({'H-optostim','V-optostim'},'location','northwest');
hTitle=title('\beta distributions for optostimulation conditions','FontWeight','Normal');

upFontSize(16,0.01);
hLegend.FontSize=14;
text((median(h_opto)+median(v_opto))/2 - .025, 5.85, '***','FontSize',20) %-0.025 to account for length of text string

%text((median(h_opto)+median(bl_opto))/2, 7.8, 'n.s.','FontSize',20)

%text((median(v_opto)+median(bl_opto))/2, 7.4, '***','FontSize',20)

if saveFlag
    saveFigure(mainPath,'Chip/Meta/summary/betaHistogram','');
end
end