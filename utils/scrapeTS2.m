function dataTable = scrapeTS2(dataFolder, monkeyName, dateStart, dateEnd)
    %% Description:
    % See original detailed description in your script.
    
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
                
                % Define the file pattern using wildcards
                filePattern = fullfile(dateFolderPath, runFolderName, ['M28D*R*TS.mat']);
                matchingFiles = dir(filePattern);
                
                for j = 1:length(matchingFiles)
                    filePath = fullfile(matchingFiles(j).folder, matchingFiles(j).name);
                    
                    % Check if the file exists
                    if isfile(filePath)
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
                        prf = 0; % New column for PRF
                        
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
                                    
                                    % Check for the PRF condition
                                    if all(cellfun(@isempty, projImages) | ...
                                            cellfun(@(x) contains(x, 'PK\LumSeries\Dot0125P00000L'), projImages))
                                        prf = 1;
                                        
                                        % Check for *StabDFFTAmp*.mat in the folder
                                        stabFiles = dir(fullfile(dateFolderPath, runFolderName, '*StabDFFTAmp*.mat'));
                                        if ~isempty(stabFiles)
                                            prf = prf + 1; % Increment PRF to indicate 1+
                                        end
                                    end
                                    
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
                        results = [results; {dateStr, runNum, protocolName, gaborContrasts, gaborOrientation, optoOnly, blOnly, bloptoCombined, intgMat, binned, prf}];
                    end
                end
            end
        end
    end
    
    % Convert results to table format
    dataTable = cell2table(results, 'VariableNames', {'Date', 'RunNumber', 'ProtocolName', 'GaborContrasts', 'GaborOrientation', 'OptoOnly', 'BlOnly', 'BlOptoCombined', 'IntgMat', 'Binned', 'PRF'});
    
    % Sort the table by Date and RunNumber in ascending order
    dataTable = sortrows(dataTable, {'Date', 'RunNumber'});
    
    % Display the table as output
    disp(dataTable);
end
