function plotOptoRF
 %% Plot
            % sort according to powerr
            [~,sortOrder]=sort(plotImg.temporalON);
            plotImg.temporalON=plotImg.temporalON(sortOrder);
            plotImg.fullROI=plotImg.fullROI(:,:,:,sortOrder);
            plotImg.recruitROI=plotImg.recruitROI(:,:,:,sortOrder);
            plotImg.fullROIcolumnar=plotImg.fullROIcolumnar(:,:,:,sortOrder);
            plotImg.recruitROIcolumnar=plotImg.recruitROIcolumnar(:,:,:,sortOrder);

            figure('name','Raw full ROI')
            nCols=size(plotImg.fullROI,4); nRows=size(plotImg.fullROI,3) + 1;
            gap=.02;marginV=.025;marginH=.01;
            [hAx,~]=tight_subplot(nRows,nCols,[gap+.02 gap], [marginV+.1 marginV+.2], [marginH+.05 marginH]);
            for row=1:nRows
                for col=1:nCols
                    plotNo=col + (row-1)*nCols;
                    axes(hAx(plotNo))

                    if row<=2
                        imgsc(plotImg.fullROI(:,:,row,col));
                        caxis([0 3.5] * 10^-1)
                    elseif row==3
                        imgsc(plotImg.fullROI(:,:,2,col)-plotImg.fullROI(:,:,1,col));
                        caxis([-8 8] * 10^-2);
                    end
                    if row==1
                        title(sprintf('%.0f %%',plotImg.temporalON(col)))
                    end
                    addPix2MM(plotImg.fullROI(:,:,1,col),plotNo,nRows,nCols);
                    colorbar; 
                end
            
                upFontSize(28,.008); axis square
            end
end