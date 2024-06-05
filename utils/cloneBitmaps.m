function bmp=cloneBitmaps(currentBlockStruct, cloneLoadFlag)
%COPYUNIQUEBMPFILES Loads a file, extracts unique BMP filenames, and copies them.
%   sourceFilePath: Path to the file containing the TS variable.
%   directoryPath: Directory where the BMP files are located.
%   sessionPath: Destination directory for the copied BMP files.
bmp=nan(1080,1920,2);
bmpPath=['T:\'];
savePath=currentBlockStruct.ROITC;

% Load the file containing the TS variable
loadedData = load(currentBlockStruct.TS);
TS = loadedData.TS;

% Extract the unique .bmp filenames from TS.Header.Conditions.ProjImg
bmpFilenamesRaw = unique(TS.Header.Conditions.ProjImg(cellfun('isempty', strfind(TS.Header.Conditions.ProjImg, 'Dot'))));
bmpFilenamesRaw=bmpFilenamesRaw(end-1:end);
bmpFilenamesRaw=strrep(bmpFilenamesRaw,'\','/');
bmpFilenames = findMatchingFiles(bmpPath,bmpFilenamesRaw);
% Ensure the sessionPath directory exists
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

switch cloneLoadFlag
    case {'clone'}
        % Look for these filenames in directoryPath and copy them to sessionPath
        for i = 1:length(bmpFilenames)
            srcFile = bmpFilenames{i};
            if exist(srcFile, 'file')
                [~, name, ext] = fileparts(bmpFilenames{i});
                destFile = fullfile(savePath, [name ext]);
                copyfile(srcFile, destFile);
                fprintf('Cloning... Copied %s to %s\n',[name ext],savePath);
            else
                fprintf('Cloning... File %s not found in %s\n', bmpFilenames{i}, bmpPath);
            end
        end
    case {'load'}
        % Look for these filenames in directoryPath and copy them to sessionPath
        for i = 1:length(bmpFilenames)
            [~, name, ext] = fileparts(bmpFilenames{i});
            destFile = fullfile(savePath, [name ext]);
            if exist(destFile, 'file')
                bmp(:,:,i)=double(imread(destFile));
                fprintf('Loading... Loaded %s\n', destFile);
            else
                bmp(:,:,i)=nan(1080,1920);
                fprintf('Loading... File %s not found\n', destFile);
            end
        end
end
end
%{
function matchingFiles = findMatchingFiles(bmpPath,filePaths)
    % Initialize output
    matchingFiles = {};
    
    % Loop over each file path
    for i = 1:length(filePaths)
        % Extract directory and filename from the file path
        [fileDir, fileName, fileExt] = fileparts([bmpPath filePaths{i}]);
        
        % Create a pattern to match the filename without the potential Cxx suffix
        filePattern = [fileName(1:end-4), '*', fileExt];  % e.g., O09000HE0032G040S00001*
        
        % Get list of matching files in the directory
        files = dir(fullfile(fileDir, filePattern));
        
        % Filter out non-bmp files
        files = files(endsWith({files.name}, '.bmp'));
        
        % Store the full paths of matching files
        for j = 1:length(files)
            matchingFiles{end+1} = fullfile(fileDir, files(j).name);
        end
    end
end
%}

function matchingFiles = findMatchingFiles(bmpPath, filePaths)
    % Initialize output
    matchingFiles = {};
    
    % Loop over each file path
    for i = 1:length(filePaths)
        % Extract directory and filename from the file path
        [fileDir, fileName, fileExt] = fileparts(fullfile(bmpPath, filePaths{i}));
        
        % Create a pattern to match the filename without the potential Cxx suffix
        filePattern = [fileName(1:end-4), '*', fileExt];  % e.g., O09000HE0032G040S00001*
        
        % Get list of matching files in the directory
        files = dir(fullfile(fileDir, filePattern));
        
        % Filter out non-bmp files
        files = files(endsWith({files.name}, '.bmp'));
        
        % Check for direct match first
        directMatch = [];
        for j = 1:length(files)
            if strcmp(files(j).name, [fileName, fileExt])
                directMatch = files(j).name;
                break;
            end
        end
        
        % If a direct match is found, add it to the matchingFiles
        if ~isempty(directMatch)
            matchingFiles{end+1} = fullfile(fileDir, directMatch);
        else
            % Find the closest match based on the last few characters
            bestMatch = '';
            bestMatchDist = Inf;
            
            for j = 1:length(files)
                currentFile = files(j).name;
                
                % Ensure both strings are of the same length for comparison
                maxLength = max(length(fileName), length(currentFile));
                paddedFileName = pad(fileName, maxLength, 'right');
                paddedCurrentFile = pad(currentFile, maxLength, 'right');
                
                % Calculate the difference in the last few characters
                dist = sum(paddedFileName ~= paddedCurrentFile);
                
                % Update the best match if the current file is a closer match
                if dist < bestMatchDist
                    bestMatch = currentFile;
                    bestMatchDist = dist;
                end
            end
            
            % Save the best match to the output cell array
            if ~isempty(bestMatch)
                matchingFiles{end+1} = fullfile(fileDir, bestMatch);
            end
        end
    end
end
