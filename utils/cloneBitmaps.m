function bmp=cloneBitmaps(currentBlockStruct, cloneLoadFlag)
% Check if TS file exists first
if ~isfile(currentBlockStruct.TS)
    % If TS file doesn't exist, return NaN matrix
    bmp = nan(1080, 1920, 2);
    return;
end

% Rest of the original function...
bmp = nan(1080, 1920, 2);
bmpPath = currentBlockStruct.mainPath;
savePath = currentBlockStruct.ROITC;

% Load the file containing the TS variable
loadedData = load(currentBlockStruct.TS);
TS = loadedData.TS;

% Extract the unique .bmp filenames from TS.Header.Conditions.ProjImg
bmpFilenamesRaw = unique(TS.Header.Conditions.ProjImg(cellfun('isempty', strfind(TS.Header.Conditions.ProjImg, 'Dot'))));
bmpFilenamesRaw=bmpFilenamesRaw(end-1:end);
bmpFilenamesRaw=strrep(bmpFilenamesRaw,'\','/');

% Ensure the sessionPath directory exists
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

switch cloneLoadFlag
    case {'clone'}
        resultArray = cellfun(@(x) extractAfter(x, 'PK'), bmpFilenamesRaw, 'UniformOutput', false);
        bmpFilenames = findMatchingFiles([bmpPath 'users/PK/colStimPipeline'],resultArray,currentBlockStruct.date);
        % Look for these filenames in directoryPath and copy them to sessionPath
        for i = 1:length(bmpFilenames)
            srcFile = bmpFilenames{i};
            if exist(srcFile, 'file')
                [~, name, ext] = fileparts(bmpFilenames{i});
                destFile = fullfile(savePath, [name ext]);
                copyfile(srcFile, destFile);
                fprintf('Cloning... Copied %s to %s\n',[name ext],savePath);
                bmp(:,:,i)=double(imread(destFile));
            else
                fprintf('Cloning... File %s not found in %s\n', bmpFilenames{i}, bmpPath);
                bmp(:,:,i)=nan(1080,1920);
            end
        end
    case {'load'}
        for i=1:length(bmpFilenamesRaw)
            [filepath,filename,fileext]=fileparts(bmpFilenamesRaw{i});

            fileNames{i}=[filename fileext];
        end
        
        bmpFilenames = findMatchingFiles(currentBlockStruct.ROITC,fileNames, currentBlockStruct.date);

        % Look for these filenames in directoryPath and copy them to sessionPath
        for i = 1:length(bmpFilenames)
            [~, name, ext] = fileparts(bmpFilenames{i});
            destFile = bmpFilenames{i}; % fullfile(savePath, [name ext])
            if exist(destFile, 'file')
                bmp(:,:,i)=double(imread(destFile));
                %fprintf('Loading... Loaded %s\n', destFile);
            else
                bmp(:,:,i)=nan(1080,1920);
                %fprintf('Loading... File %s not found\n', destFile);
            end
        end
end
end

function matchingFiles = findMatchingFiles(bmpPath, filePaths, currentDate)
    % Initialize output
    matchingFiles = {};

    % Loop over each file path
    for i = 1:length(filePaths)
        % Extract directory and filename from the file path
        [fileDirData, fileName, fileExt] = fileparts(fullfile(bmpPath, filePaths{i}));

        % Initialize a flag to indicate if a match has been found
        foundMatch = false;

        % Define the directories to search
        searchDirs = {fileDirData,['V:/PK/ColSeries/' currentDate]}; % Search in either the server Y:\ or local experiment computer T:\

        % Iterate over search directories
        for k = 1:length(searchDirs)
            % Create a pattern to match the filename without the potential Cxx suffix
            filePattern = [fileName(1:end-8), '*', fileExt];  % e.g., O09000HE0032G040S00001*

            % Get list of matching files in the current directory
            files = dir(fullfile(searchDirs{k}, filePattern));

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
                matchingFiles{end+1} = fullfile(searchDirs{k}, directMatch);
                foundMatch = true;
                break; % Stop searching once a match is found
            end
        end

        % If no direct match was found in any directory, add logic for closest match if needed
        if ~foundMatch
            % Repeat search logic for closest match (if applicable)
            % Similar code block can be added here if no matches are found
        end
    end
end

%{
function matchingFiles = findMatchingFiles(bmpPath, filePaths)
    % Initialize output
    matchingFiles = {};
    
    % Loop over each file path
    for i = 1:length(filePaths)
        % Extract directory and filename from the file path
        [fileDir, fileName, fileExt] = fileparts(fullfile(bmpPath, filePaths{i}));
        
        % Create a pattern to match the filename without the potential Cxx suffix
        filePattern = [fileName, '*', fileExt];  % e.g., O09000HE0032G040S00001*
        
        % Get list of matching files in the directory
        files = dir(fullfile([bmpPath fileDir(end-18:end)], filePattern));
        
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
%}
