%{
datastruct.DCseries={};
datastruct.SIRFseries={};
datastruct.PRFseries={};
%}



% Define the base directory containing the data
baseDir = 'Y:\Chip';

% Construct a search pattern for files matching the format
searchPattern = fullfile(baseDir, 'Chip2024*\', 'run*\', 'M*TS.mat');
%load('Y:\Chip\forGiac\Chip20220422\run2\M28D20220422R2TS.mat')
%TS.Header.Conditions.ProjImg'

% Find all files matching the pattern
allFiles = dir(searchPattern);

% Initialize an empty structure to store filenames
dataStruct = struct('LumSeries', {}, 'DCseries', {}, 'SFseries', {});

for iFile = 1:length(allFiles)
  % Get current file information
  fileInfo = allFiles(iFile);
  
  % Construct the full file path
  filePath = fullfile(fileInfo.folder, fileInfo.name);
  
  % Load the data
  try
    loadedData = load(filePath);
  catch ME
    % Handle potential errors during loading
    warning(['Error loading file: ' filePath]);
    continue;
  end
  
  % Check for required fields in the loaded data
  if ~isfield(loadedData, 'TS') || ~isfield(loadedData.TS, 'Header') || ...
      ~isfield(loadedData.TS.Header, 'Conditions') || ~isfield(loadedData.TS.Header.Conditions, 'ProjImg')
    continue;
  end
  
  % Get the ProjImg value
  projImg = loadedData.TS.Header.Conditions.ProjImg;


    % Check if any element in projImg contains the desired substrings
    if any(cellfun(@(x) any(contains(x, {'LumSeries', 'DCseries', 'SFseries', 'Checkboard'})), projImg))
        if sum(contains(projImg,'LumSeries'))>0
            datastruct.PRFseries(end+1,1)={filePath};
            disp(['File: ' filePath ' belongs to PRF series (' projImg{2} ')']);
        end
        
        if sum(contains(projImg,'DCseries'))>0
            datastruct.DCseries(end+1,1)={filePath};
            disp(['File: ' filePath ' belongs to DC series (' projImg{2} ')']);
        end
        
        if sum(contains(projImg,'SFseries'))>0
            datastruct.SIRFseries(end+1,1)={filePath};
            disp(['File: ' filePath ' belongs to SF series (' projImg{2} ')']);
        end

    end

end

%Display the contents of the data structure (filenames)
disp('Files in LumSeries:')
disp(datastruct.PRFseries)
disp('Files in DCseries:')
disp(datastruct.DCseries)
disp('Files in SFseries:')
disp(datastruct.SIRFseries)

save(['Y:/Chip/Meta/summary/ChipPatternedOptostimPaper.mat'],'datastruct')