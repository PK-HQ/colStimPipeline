
function generateOptostimTable(chamberWanted, bitmapData)

% Extract entries for given chamber
[mainPath, datastruct]=setupEnv('users/PK/colStimPipeline/exptListBiasingFull.m');
nColumnsWanted=[]; 
chamberEntries = organizeBlocks(datastruct, chamberWanted, nColumnsWanted);

% Extract datastruct entries for given chamber
exptList=datastruct(chamberEntries);

% Add fields from bitmapData to exptList
fieldsToAdd = {'ledpower_mW', ...
               'areaOrtMask', 'areaGaussMask', 'areaPixelsON',... 
               'pixelsON', 'timeONPercent', 'adjustedSPD_uW'};
exptList = addFieldsToStruct(exptList, bitmapData, fieldsToAdd);

% Convert exptBiasingListFull datastruct to table
dataTable=struct2table(exptList);
noTableHeaders=size(dataTable,2);
sortIdx=[1:4 19 21 22 5 6 7  23 9 8 11:15 12 13 15 18 16 17 20 24:noTableHeaders]; % sort for excel table by imaging and optostim
dataTable=dataTable(:,sortIdx);

% Add new column with chamber name in front of table
numRows = height(dataTable);  % Gets the number of rows in the table T
chamberColumn = repmat({chamberWanted}, numRows, 1);  % Create a cell array with the repeated string
newColumnTable = table(chamberColumn, 'VariableNames', {'Chamber'});
dataTable = [newColumnTable dataTable];

excelFilename=[mainPath 'users\PK\colStimPipeline\docs\optostimTable'];
writetable(dataTable,[excelFilename '06202024.xlsx'],'Sheet',chamberWanted)

end

function originalStruct = addFieldsToStruct(originalStruct, addedStruct, fieldsToAdd)
    % Ensure the number of elements in exptList matches the columns in each field of bitmapData
    for k = 1:length(fieldsToAdd)
        fieldName = fieldsToAdd{k};
        fieldData=squeeze(addedStruct.(fieldName))';
        for row=1:size(fieldData,1)
            originalStruct(row).(fieldName)=fieldData(row,:);
        end
    end
end
