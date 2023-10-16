function addPix2MM(minXpx,maxXpx,minYpx,maxYpx,subpltIdx,nRows,nCols)
% Imaging window in pixels and milimeters
mmTotal=8.22;
pxTotal=512;
% Conversions
MMperPX=mmTotal/pxTotal;
PXperMM=pxTotal/mmTotal;
% X-Y limits in pixels
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

if subpltIdx==selectedSubplot
    
    %xlabel('X (mm)');
    %ylabel('Y (mm)');

    % add xticks
    xticks(minXpx-minXpx:PXperMM:maxXpx-minXpx)
    nXMM=maxXMM-minXMM;
    % add xticklabels
    xlabels = string([0:1:nXMM]); % define xticklabels
    xlabels(2:2:end) = NaN; % remove every other one
    ax.XAxis.TickLabels = xlabels; % set
    % set xlim
    xlim([0 maxXpx-minXpx])

    %same for y
    yticks(minYpx-minYpx:PXperMM:maxYpx-minYpx)
    nYMM=maxYMM-minYMM;
    ylabels = string([0:1:nYMM]); % extract
    ylabels(2:2:end) = NaN; % remove every other one
    ax.YAxis.TickLabels = ylabels; % set
    ylim([0 maxYpx-minYpx])
else
    % add xticks
    xticks(minXpx-minXpx:PXperMM:maxXpx-minXpx)
    nXMM=maxXMM-minXMM;
    % add xticklabels
    xlabels = string([0:1:nXMM]); % define xticklabels
    xlabels(:) = NaN; % remove every other one
    ax.XAxis.TickLabels = xlabels; % set
    % set xlim
    xlim([0 maxXpx-minXpx])

    %same for y
    yticks(minYpx-minYpx:PXperMM:maxYpx-minYpx)
    nYMM=maxYMM-minYMM;
    ylabels = string([0:1:nYMM]); % extract
    ylabels(:) = NaN; % remove every other one
    ax.YAxis.TickLabels = ylabels; % set
    ylim([0 maxYpx-minYpx])
end
end