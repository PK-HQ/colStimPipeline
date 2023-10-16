function pubfig(varargin)

if length(varargin)==1;
  figHandle=varargin{1};
  fontSize=12;
  lineWeight=3;
  tickLength=.5;
else
  figHandle=varargin{1};
  fontSize=varargin{2};
  lineWeight=varargin{3};
  tickLength=varargin{4};
end

%global
%set(groot,'defaultLineLineWidth',lineWeight)
%set(0,'DefaultAxesTitleFontWeight','normal');

% Set title to non-bold, and 1.3x font size of other text size
set(0,'DefaultAxesTitleFontWeight','normal');
%set(0,'DefaultAxesTitleFontSize',1.5);

%change font sizes of all elements
set(findall(gcf,'-property','FontSize'),'FontSize',fontSize)
%retain only axes
box off

%set linewidth of axes and plot
set(gca,'LineWidth',lineWeight); %ax
if isprop(figHandle,'LineWidth')==1
    set(figHandle,'LineWidth',lineWeight); %plot
end
%ticks face outwards
set(gca,'TickDir','out');

%ticks longer
ax = gca;
ax.TickLength = [tickLength,tickLength];
%white background
set(gcf,'color','w')

end