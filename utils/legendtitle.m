function legendtitle(legendTitle)
    hLeg = findobj(gcf, 'Type', 'Legend');
    hLegTitle = get(hLeg(1),'Title');
    fontsize = get(hLeg(1),'FontSize');
    set(hLegTitle,'String',legendTitle,'FontSize',fontsize);
end
