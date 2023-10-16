function [dotProduct,cosAngle,corrValue,intensityPrct]=calculateSimilarity(img,imgReference)
%dot product
dotProduct=sum(img.*imgReference,'all');
cosAngle=dotProduct/(norm(img)*norm(imgReference));

%correlation r-value
corrMatrix=corrcoef(img,imgReference,'rows','complete');
corrValue=corrMatrix(2);

% Crude intensity
%{
intensityDiff=nansum(img,'all')-nansum(imgReference(:),'all');
intensityPrct=100 * (intensityDiff./nansum(imgReference(:),'all'));
%}
 
% Spatial intensity
intensityDiff=img-imgReference;
intensityPrct=intensityDiff*100./imgReference;
intensityPrct(isinf(intensityPrct))=NaN;
intensityPrct=nanmean(intensityPrct(:));

end