function histEq=histEqualization(normVERpca,imgDims,gridSize)
for ortNo=1:size(normVERpca,3)
  histEq(:,:,ortNo)=adapthisteq(normVERpca(:,:,ortNo),'NumTiles',imgDims./gridSize,'Range','full');
end
end