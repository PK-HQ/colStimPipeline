function gammaCorrection(histEqCol,gammaCorr)
for ortNo=1:size(histEqCol,3)
  gammaCol(:,:,ortNo)=imadjust(histEqCol(:,:,ortNo),[],[],gammaCorr);
end