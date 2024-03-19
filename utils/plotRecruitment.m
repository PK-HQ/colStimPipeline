function plotRecruitment(plotImg, roiType, nRows, nCols, caxLimit, titleStr)
gap = .02; marginV = .025; marginH = .01;
[hAx, ~] = tight_subplot(nRows, nCols, [gap+.02 gap], [marginV+.1 marginV+.15], [marginH+.05 marginH+.025]);
rowLabels={'0\circ','90\circ','90-0\circ'};

for row = 1:nRows
    for col = 1:nCols
        plotNo = col + (row - 1) * nCols;
        axes(hAx(plotNo));

        if row <= 2
            imgsc(plotImg.(roiType)(:,:,row,col));
            cax(caxLimit, plotImg.(roiType)(:,:,1:2,:));
        elseif row == 3
            imgsc(plotImg.(roiType)(:,:,2,col) - plotImg.(roiType)(:,:,1,col));
            cax(caxLimit, plotImg.(roiType)(:,:,2,:) - plotImg.(roiType)(:,:,1,:));
        end

        % Titles and labels
        addPix2MM(plotImg.(roiType)(:,:,1,col), plotNo, nRows, nCols);
        if row == 1
            title(sprintf('%.0f %%', plotImg.temporalON(col)));
        end

        addVerticalLabel(hAx, row, col, nCols, rowLabels{row});

        %axes
        colorbar;
    end
end
[ax,h]=suplabel(titleStr,'t',[0.08 0.08 .87 .86]);
set(h,'FontWeight','normal');
upFontSize(28, .008); axis square;
end