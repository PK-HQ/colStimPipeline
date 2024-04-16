function resultTable = convertStruct2Table(structData, fieldName)
    % Load the MAT file structName
    %structName = load(matFile);
    % Assume the struct name is bitmapData and access it
    %s%tructData = structData.bitmapData;

    % Prepare to capture the structName in cell arrays
    numBlocks = size(structData, 3);
    tableData = cell(numBlocks, 1);

    % Loop through each block
    for i = 1:numBlocks
        % Access the field for both bitmaps and concatenate them into a string
        % Check if we need to squeeze any dimension
        value1 = squeeze(structData(:, i).(fieldName));
        value2 = squeeze(structData(:, i).(fieldName));
        
        % Ensure the structName is not empty and convert to string
        if isempty(value1)
            val1Str = '';
        else
            val1Str = num2str(value1);
        end
        
        if isempty(value2)
            val2Str = '';
        else
            val2Str = num2str(value2);
        end

        % Combine both strings separated by a comma
        tableData{i} = [val1Str, ', ', val2Str];
    end

    % Create a table from the cell array
    resultTable = table(tableData, 'VariableNames', {fieldName});
end
