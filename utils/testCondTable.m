condTableSize=size(TS.Header.ConditionTable,1);
for i=8:condTableSize
    inputStr = TS.Header.ConditionTable(8);
    numSpaces = 10; % Number of additional spaces to add
    formattedStr(i) = addSpaces(inputStr, numSpaces);
end

