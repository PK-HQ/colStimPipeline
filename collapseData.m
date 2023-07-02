function [uniqueOrts, correctPercent, errorPercent]=collapseData(ortData, correctData)
uniqueOrts=sort(unique(ortData),'descend'); %smallest on top
for ortNo=1:numel(uniqueOrts)
    ort=uniqueOrts(ortNo);
    outcomes=correctData(find(ortData==ort));
    outcomes=outcomes;
    correctPercent(ortNo)=sum(outcomes) * 100 ./numel(outcomes);
    errorPercent(ortNo)=100-correctPercent(ortNo);
end
end