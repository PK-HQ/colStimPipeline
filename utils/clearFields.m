function datastruct=clearFields(datastruct, fieldsToKeep)
% Create a new datastruct with only the fields to keep
newdatastruct = struct();
for k = 1:length(fieldsToKeep)
    fieldName = fieldsToKeep{k};
    if isfield(datastruct, fieldName)
        newdatastruct.(fieldName) = datastruct.(fieldName);
    end
end

% Replace the original datastruct with the new one
datastruct = newdatastruct;

end