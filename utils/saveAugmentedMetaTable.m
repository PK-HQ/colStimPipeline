function saveAugmentedMetaTable(MetaTable, datastruct, analysisBlockID, columnsDesired)
% Overwrite the standard MetaTable MAT/XLSX files after enrichment.

if ~isstruct(datastruct) || isempty(datastruct) || ...
        ~isfield(datastruct, 'monkey') || ~isfield(datastruct, 'chamber')
    warning('MetaTable:SaveSkipped', ...
        'Could not infer monkey/chamber; augmented MetaTable was not saved.');
    return
end

entry = datastruct(1);
if ~isempty(analysisBlockID) && analysisBlockID(1) >= 1 && ...
        numel(datastruct) >= analysisBlockID(1)
    entry = datastruct(analysisBlockID(1));
end

outDir = ['Y:/' entry.monkey '/Meta/summary'];
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
baseName = ['metaTable-' entry.chamber num2str(columnsDesired)];
save([outDir '/' baseName '.mat'], 'MetaTable');

xlsxPath = [outDir '/' baseName '.xlsx'];
excelCell = tableToExcelCell(MetaTable);
if exist(xlsxPath, 'file')
    delete(xlsxPath);
end
writecell(excelCell, xlsxPath, 'Sheet', 'metadata');
end


function excelCell = tableToExcelCell(MetaTable)
headers = MetaTable.Properties.VariableDescriptions;
if numel(headers) ~= width(MetaTable) || any(cellfun(@isempty, headers))
    headers = MetaTable.Properties.VariableNames;
end
excelCell = cell(height(MetaTable) + 1, width(MetaTable));
excelCell(1, :) = headers;

for row = 1:height(MetaTable)
    for col = 1:width(MetaTable)
        name = MetaTable.Properties.VariableNames{col};
        if iscell(MetaTable.(name))
            value = MetaTable.(name){row};
        else
            value = MetaTable.(name)(row, :);
        end
        excelCell{row + 1, col} = excelScalarize(value);
    end
end
end


function output = excelScalarize(value)
if isempty(value)
    output = '';
elseif (isnumeric(value) || islogical(value)) && isscalar(value)
    output = value;
elseif isnumeric(value) || islogical(value)
    if ndims(value) <= 2
        output = mat2str(value, 6);
    else
        output = sprintf('<%s [%s], min=%g, max=%g, mean=%g>', ...
            class(value), strtrim(sprintf('%d ', size(value))), ...
            min(value(:), [], 'omitnan'), max(value(:), [], 'omitnan'), ...
            mean(double(value(:)), 'omitnan'));
    end
elseif ischar(value)
    output = value;
elseif isstring(value)
    output = char(strjoin(cellstr(value), '; '));
else
    try
        output = char(string(value));
    catch
        output = ['<' class(value) '>'];
    end
end
end
