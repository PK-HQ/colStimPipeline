function plotCounter=plotter(dataMat,mask,selectedOrts,plotIdx,plotCounter,titleStr,nRows,nCols)
Ort=[0 90];
imgDims=size(dataMat(:,:,1));
imgX=[0 imgDims(1)];
imgY=[0 imgDims(2)]
for ortNo=1:length(selectedOrts)
  axes(hAx(plotIdx(plotCounter)));
  img=dataMat(:,:,ortNo);
  greyNaN(img.*mask) %.*ROItuningMap optional for raw
  titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
  title({titleStr,titleL1},'FontWeight','normal');
  axis square;
  addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
  plotCounter=plotCounter+1;
  
  if plotCounter==1;
    xlabel('X (mm)');
    ylabel('Y (mm)');
  end
end





;


upFontSize(14,.015)
set(h,'FontSize',16)
end
