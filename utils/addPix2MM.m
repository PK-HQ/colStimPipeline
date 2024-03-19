function addPix2MM(image,subplotIdx,nRows,nCols)
% Imaging window in pixels and milimeters
mmTotal=8.22;
pxTotal=512;
% Conversions
MMperPX=mmTotal/pxTotal;
PXperMM=pxTotal/mmTotal;
% X-Y limits in pixels
minXpx=1;minYpx=1;
maxXpx=size(image,2);maxYpx=size(image,1);
minXMM=minXpx*MMperPX;
maxXMM=maxXpx*MMperPX;
minYMM=minYpx*MMperPX;
maxYMM=maxYpx*MMperPX;
ax=gca;

upFontSize(12,.02)
box off
set(gca,'TickDir','out')

%% If bottom left plot, add labels and tick marks
selectedSubplot=((nRows-1)*nCols)+1;

if subplotIdx==selectedSubplot
    
    %xlabel('X (mm)');
    %ylabel('Y (mm)');

    % add xticks
    xticks(minXpx-minXpx:PXperMM:maxXpx-minXpx)
    nXMM=maxXMM-minXMM;
    % add xticklabels
    xticklabels = string([0:1:nXMM]); % define xticklabels
    xticklabels(2:2:end) = NaN; % remove every other one
    ax.XAxis.TickLabels = xticklabels; % set
    % set xlim
    xlim([0 maxXpx-minXpx])

    %same for y
    yticks(minYpx-minYpx:PXperMM:maxYpx-minYpx)
    nYMM=maxYMM-minYMM;
    yticklabels = string([0:1:nYMM]); % extract
    yticklabels(2:2:end) = NaN; % remove every other one
    ax.YAxis.TickLabels = yticklabels; % set
    ylim([0 maxYpx-minYpx])
    xlabel('mm')
    ylabel('mm')
    
else
    set(gca,'Yticklabel',[]) 
    set(gca,'Xticklabel',[])
    %{
    % add xticks
    xticks(minXpx-minXpx:PXperMM:maxXpx-minXpx)
    nXMM=maxXMM-minXMM;
    % add xticklabels
    xticklabels = string([0:1:nXMM]); % define xticklabels
    xticklabels(:) = NaN; % remove every other one
    ax.XAxis.TickLabels = xticklabels; % set
    % set xlim
    xlim([0 maxXpx-minXpx])

    %same for y
    yticks(minYpx-minYpx:PXperMM:maxYpx-minYpx)
    nYMM=maxYMM-minYMM;
    yticklabels = string([0:1:nYMM]); % extract
    yticklabels(:) = NaN; % remove every other one
    ax.YAxis.TickLabels = yticklabels; % set
    ylim([0 maxYpx-minYpx])
    %}
end
xtickangle(0)
end