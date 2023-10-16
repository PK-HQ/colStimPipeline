function [betaOrientation,percentVertReport]=mergeZeros(betaOrientation,percentVertReport)
zeroidx=find(betaOrientation==0);
if ~isempty(zeroidx)
    meanZeros=mean(percentVertReport(zeroidx));
    percentVertReport(zeroidx(1))=meanZeros;
    percentVertReport(zeroidx(2))=NaN;
    betaOrientation(zeroidx(2))=NaN;
    percentVertReport=percentVertReport(~isnan(percentVertReport));
    betaOrientation=betaOrientation(~isnan(betaOrientation));
else
    %no change
end
end