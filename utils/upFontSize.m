function upFontSize(size,tickLength)
%INCREASE OVERALL FONT SIZE
    set(findall(gcf,'-property','FontSize'),'FontSize',size) %increase font size
    set(gca,'linewidth',2) %increase line width
    set(gca,'TickDir','out');
    set(gca,'FontWeight','Normal');
    set(gca,'TickLength',[tickLength, tickLength])
    set(gcf,'color','w');
    box off
    %remove up and right ticks
    %removeURTicks
end