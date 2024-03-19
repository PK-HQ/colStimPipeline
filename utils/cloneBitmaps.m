function bmp=cloneBitmaps(currentBlockStruct, cloneLoadFlag)
%COPYUNIQUEBMPFILES Loads a file, extracts unique BMP filenames, and copies them.
%   sourceFilePath: Path to the file containing the TS variable.
%   directoryPath: Directory where the BMP files are located.
%   sessionPath: Destination directory for the copied BMP files.

bmpPath=['T:\'];
savePath=currentBlockStruct.ROITC;

% Load the file containing the TS variable
loadedData = load(currentBlockStruct.TS);
TS = loadedData.TS;

% Extract the unique .bmp filenames from TS.Header.Conditions.ProjImg
bmpFilenames = unique(TS.Header.Conditions.ProjImg(~cellfun('isempty',TS.Header.Conditions.ProjImg)));

% Ensure the sessionPath directory exists
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

switch cloneLoadFlag
    case {'clone'}
        % Look for these filenames in directoryPath and copy them to sessionPath
        for i = 1:length(bmpFilenames)
            srcFile = fullfile(bmpPath, bmpFilenames{i});
            if exist(srcFile, 'file')
                [~, name, ext] = fileparts(bmpFilenames{i});
                destFile = fullfile(savePath, [name ext]);
                copyfile(srcFile, destFile);
                fprintf('Copied %s to %s\n', bmpFilenames{i},savePath);
            else
                fprintf('File %s not found in %s\n', bmpFilenames{i}, bmpPath);
            end
        end
    case {'load'}
        % Look for these filenames in directoryPath and copy them to sessionPath
        for i = 1:length(bmpFilenames)
            [~, name, ext] = fileparts(bmpFilenames{i});
            destFile = fullfile(savePath, [name ext]);
            if exist(destFile, 'file')
                bmp(:,:,i)=double(imread(destFile));
                fprintf('Loaded %s\n', destFile);
            else
                fprintf('File %s not found\n', destFile);
            end
        end
end
end
