function optostimOrientations=extractOptostimOrt(TS)
%% Extracts optostim pattern orientation (0 or 90deg columns) from projector image name with regex

% Define the cell array
cellArray = TS.Header.Conditions.ProjImg(find(TS.Header.Conditions.TypeCond>0));

% Preallocate a cell array to store results
optostimOrientations = nan(size(cellArray));

% Loop over each cell in the array
for i = 1:length(cellArray)
    if ~isempty(cellArray{i})
        % Use regular expression to find the pattern
        tokens = regexp(cellArray{i}, 'O0\d{2}', 'match');
        if ~isempty(tokens)
            % Store the first match (assuming there's only one match per cell)
            optostimOrtStr=tokens{1};
            optostimOrt=str2num(optostimOrtStr(2:end));
            optostimOrientations(i)=optostimOrt;
        end
    end
end
end
