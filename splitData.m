function [ortSplit1,correctSplit1,ortSplit2,correctSplit2]=splitData(dataMatrix,nSplits, splitMode)
%split sorted data (chronological) into 2 parts (early, late)
ortAll=dataMatrix(:,2);
correctAll=dataMatrix(:,5);


switch splitMode
    case {'trial'} %by trial number, preserves time of trial presented but will end up with odd numbers of trials per cond 
        switch nSplits
            case {2}
                %early
                ortSplit1=ortAll(1:floor(size(ortAll,1)/2));
                correctSplit1=correctAll(1:floor(size(correctAll,1)/2));
                %late
                ortSplit2=ortAll(floor(size(ortAll,1)/2)+1:end);
                correctSplit2=correctAll(floor(size(correctAll,1)/2)+1:end);
            case {1}
                ortSplit1=ortAll;
                correctSplit1=correctAll;
                ortSplit2=[];
                correctSplit2=[];
        end
    case {'cond'} %by trial number, preserves n-trials per cond but destroys time of trial presented
        

end
