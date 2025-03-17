function dataTable = scrapeTS(dataFolder, monkeyName, dateStart, dateEnd)
 %% Description:
% This function loads data from specific MATLAB `.mat` files located in a nested folder structure. 
% The function accepts a base data folder, a monkey name, and a date range as inputs. It searches
% within the specified date range for files matching a specific pattern. If a file contains 
% certain variables, their values are extracted and added to a results table.
% 
% The function dynamically locates all run folders within each date folder and processes only
% existing files, making it efficient and adaptable to variable folder contents.
%
% Parameters:
% - dataFolder (string): The base directory containing the data files (e.g., 'Y:\').
% - monkeyName (string): The name of the monkey, used as a folder identifier (e.g., 'Chip\').
% - dateStart (string): The start date in 'YYYYMMDD' format. Only files from this date onward 
%   will be processed.
% - dateEnd (string): The end date in 'YYYYMMDD' format. Only files up to this date will be 
%   processed.
%
% Returns:
% - dataTable (table): A MATLAB table containing the following columns:
%   - Date (string): The date (YYYYMMDD) associated with each processed file.
%   - RunNumber (integer): The run number extracted from each run folder.
%   - ProtocolName (string): The protocol name from `TS.Header.ProtocolName`, or an error message
%     if the file could not be loaded.
%   - GaborContrasts (array): Unique values from `TS.Header.Conditions.StimCon` if available.
%   - GaborOrientation (array): Unique values from `TS.Header.Conditions.GaborOrt` if available.
%   - OptoOnly (logical): Set to 1 if only opto condition images ('O000' and 'O090') are found 
%     in `TS.Header.Conditions.ProjImg` without 'Dot'.
%   - BlOnly (logical): Set to 1 if only baseline condition images ('Dot') are found in 
%     `TS.Header.Conditions.ProjImg`.
%   - BlOptoCombined (logical): Set to 1 if both baseline and opto conditions are found 
%     ('Dot', 'O000', and 'O090') in `TS.Header.Conditions.ProjImg`.
%
% Behavior:
% - For each date in the specified range, the function checks if a folder named `ChipYYYYMMDD` 
%   exists in `[dataFolder monkeyName]`. 
% - It then locates all run folders within the date folder, identified by folders named 'run#'.
% - For each run folder, the function attempts to load a `.mat` file following the pattern 
%   `M28YYYYMMDDR#TS.mat`, where YYYYMMDD is the date, and R# is the run number.
% - If the file exists, it is loaded, and the function checks for the `TS` variable and its fields:
%     - `TS.Header.ProtocolName`: Extracted as the protocol name.
%     - `TS.Header.Conditions.StimCon`: If available, unique contrast values are stored in the 
%       `GaborContrasts` column.
%     - `TS.Header.Conditions.GaborOrt`: If available, unique orientation values are stored in 
%       the `GaborOrientation` column.
%     - `TS.Header.Conditions.ProjImg`: The function searches for specific strings ('Dot', 
%       'O000', and 'O090') and categorizes the conditions in mutually exclusive columns:
%         - `OptoOnly`: Only 'O000' and 'O090' found.
%         - `BlOnly`: Only 'Dot' found.
%         - `BlOptoCombined`: All 'Dot', 'O000', and 'O090' found.
% - If the file cannot be loaded (e.g., due to corruption), the error message is stored in 
%   the `ProtocolName` column.
%
% Example:
% dataFolder = 'Y:\';
% monkeyName = 'Chip\';
% dateStart = '20230101';
% dateEnd = '20230110';
% dataTable = scrapeTS(dataFolder, monkeyName, dateStart, dateEnd);
% e.g. dataTable = scrapeTS('Y:\', 'Chip', '20221028', '20230209'), dataTable = scrapeTS('Y:\', 'Chip', '20230721', '20240202')
% This will search for data files in the specified date range for the monkey named 'Chip' and 
% output a table of protocol information, contrast values, orientation values, and conditions
% based on each valid file found.

%% Function
    % Initialize an empty array to store the results
    results = [];
    
    % Convert date range to datetime objects
    startDate = datetime(dateStart, 'InputFormat', 'yyyyMMdd');
    endDate = datetime(dateEnd, 'InputFormat', 'yyyyMMdd');
    
    % Iterate over each date within the range
    for currentDate = startDate:endDate
        % Format the current date as YYYYMMDD
        dateStr = datestr(currentDate, 'yyyymmdd');
        
        % Construct the date folder path
        dateFolderPath = fullfile(dataFolder, monkeyName, ['Chip' dateStr]);
        
        % Check if the date folder exists
        if isfolder(dateFolderPath)
            % Get a list of run folders within the date folder
            runFolders = dir(fullfile(dateFolderPath, 'run*'));
            
            % Iterate over each run folder found
            for i = 1:length(runFolders)
                runFolderName = runFolders(i).name;
                % Extract run number from the folder name (e.g., "run3" -> 3)
                runNum = sscanf(runFolderName, 'run%d');
                
                % Construct the file path for the TS file in the run folder
                %filePath = fullfile(dateFolderPath, runFolderName, ['M28D' dateStr 'R' num2str(runNum) 'TS.mat']);
                % Define the file pattern using wildcards
                filePattern = fullfile(dateFolderPath, runFolderName, ['M28D*R*TS.mat']);

                % Find files matching the pattern
                matchingFiles = dir(filePattern);

                % Check if any files are found
                if ~isempty(matchingFiles)
                    for i=1:length(matchingFiles)
                        % Load the first matching file (or handle multiple files as needed)
                        filePaths{i} = fullfile(matchingFiles(i).folder, matchingFiles(i).name);
                    end
                end
                for i=1:length(filePaths)
                    filePath=filePaths{i};
                    % Check if the file exists
                    if isfile(filePath)
                        % Print the file being loaded
                        fprintf('Loading file %s...\n', filePath);

                        % Initialize variables for the output table
                        protocolName = '';
                        gaborContrasts = [];
                        gaborOrientation = [];
                        optoOnly = 0;
                        blOnly = 0;
                        bloptoCombined = 0;
                        intgMat = 0;
                        binned = 0;

                        try
                            % Attempt to load the file and check if TS variable exists
                            loadedData = load(filePath, 'TS');

                            if isfield(loadedData, 'TS')
                                % Extract ProtocolName
                                if isfield(loadedData.TS, 'Header') && isfield(loadedData.TS.Header, 'ProtocolName')
                                    protocolName = loadedData.TS.Header.ProtocolName;
                                end

                                % Check for gaborContrasts (StimCon)
                                if isfield(loadedData.TS.Header, 'Conditions') && isfield(loadedData.TS.Header.Conditions, 'StimCon')
                                    gaborContrasts = unique(loadedData.TS.Header.Conditions.StimCon);
                                end

                                % Check for gaborOrientation (GaborOrt)
                                if isfield(loadedData.TS.Header.Conditions, 'GaborOrt')
                                    gaborOrientation = unique(loadedData.TS.Header.Conditions.GaborOrt);
                                end

                                % Check for ProjImg and categorize opto conditions
                                if isfield(loadedData.TS.Header.Conditions, 'ProjImg')
                                    projImages = loadedData.TS.Header.Conditions.ProjImg;
                                    containsDot = any(contains(projImages, 'Dot'));
                                    containsO000 = any(contains(projImages, 'O000'));
                                    containsO090 = any(contains(projImages, 'O090'));

                                    % Set the mutually exclusive columns
                                    if containsO000 && containsO090 && ~containsDot
                                        optoOnly = 1;
                                    elseif containsDot && ~containsO000 && ~containsO090
                                        blOnly = 1;
                                    elseif containsDot && containsO000 && containsO090
                                        bloptoCombined = 1;
                                    end
                                end
                            else
                                % If TS is missing, specify the absence
                                protocolName = 'TS variable missing';
                            end
                        catch ME
                            % If an error occurs, store the error message in ProtocolName
                            protocolName = ME.message;
                        end

                        % Check for presence of *Intg*.mat file in the run folder
                        intgFiles = dir(fullfile(dateFolderPath, runFolderName, '*Intg*.mat'));
                        if ~isempty(intgFiles)
                            intgMat = 1;
                        end

                        % Check for presence of *Bin*.mat file in the run folder
                        binFiles = dir(fullfile(dateFolderPath, runFolderName, '*Bin*.mat'));
                        if ~isempty(binFiles)
                            binned = 1;
                        end

                        % Append the data to results array
                        results = [results; {dateStr, runNum, protocolName, gaborContrasts, gaborOrientation, optoOnly, blOnly, bloptoCombined, intgMat, binned}];
                    end
                end
            end
        end
    end
    
    % Convert results to table format
    dataTable = cell2table(results, 'VariableNames', {'Date', 'RunNumber', 'ProtocolName', 'GaborContrasts', 'GaborOrientation', 'OptoOnly', 'BlOnly', 'BlOptoCombined', 'IntgMat', 'Binned'});
    
    % Sort the table by Date and RunNumber in ascending order
    dataTable = sortrows(dataTable, {'Date', 'RunNumber'});
    
    % Add OptoFixOnly column based on conditions
    dataTable.OptoFixOnly = strcmp(dataTable.ProtocolName, 'Flashing Grating + Projector') & ...
                            cellfun(@isempty, dataTable.GaborContrasts) & ...
                            dataTable.OptoOnly == 1;
    
    % Display the table as output
    disp(dataTable);
end
