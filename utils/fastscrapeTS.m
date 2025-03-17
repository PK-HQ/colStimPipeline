function dataTable = fastscrapeTS(dataFolder, monkeyName, dateStart, dateEnd)
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
                filePatternTbl = fullfile(dateFolderPath, runFolderName, ['M28D*R*Tables.txt']);
                matchingFilesTbl = dir(filePatternTbl);
                filePatternFFT = fullfile(dateFolderPath, runFolderName, ['M28D*R*StabDFFTAmp*.mat']);
                matchingFilesFFT = dir(filePatternFFT);

                for j = 1:length(matchingFilesTbl)
                    prf = 0;
                    if ~isempty(matchingFilesTbl)
                        filePathTbl = fullfile(matchingFilesTbl(j).folder, matchingFilesTbl(j).name);
                    end
                    if ~isempty(matchingFilesFFT)
                        filePathFFTexists = 1;
                    else
                        filePathFFTexists = 0;
                    end

                    % Check if the file exists
                    if isfile(filePathTbl)
                        fprintf('Loading file %s...\n', filePathTbl);
                    
                        % Read the file as text
                        fid = fopen(filePathTbl, 'r');
                        rawData = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
                        fclose(fid);
                        rawLines = rawData{1};
                    
                        % Concatenate all lines into a single string
                        concatenatedLines = strjoin(rawLines, '\n');
                    
                        % Define the search string with escaped backslashes
                        stringToFind = 'PK\\LumSeries\\Dot0125P00000L100';
                    
                        % Search for the string in the concatenated text
                        regexpMatch = ~isempty(regexp(concatenatedLines, stringToFind, 'once'));
                    
                        % Determine PRF value
                        if regexpMatch && filePathFFTexists
                            prf = 2;
                        elseif regexpMatch
                            prf = 1;
                        else
                            prf = 0;
                        end
                    
                        % Append the data to results array
                        results = [results; {dateStr, runNum, prf}];
                    end
                end
            end
        end
    end
    
    % Convert results to table format
    dataTable = cell2table(results, 'VariableNames', {'Date', 'RunNumber', 'PRF'});
    
    % Sort the table by Date and RunNumber in ascending order
    dataTable = sortrows(dataTable, {'Date', 'RunNumber'});
    
    % Display the table as output
    disp(dataTable);
end
