function [analysisBlockID] = organizeBlocks(varargin)
if length(varargin)==2
    datastruct=varargin{1};
    chamberWanted=varargin{2};
    nColumnsWanted=[];
    monkeyName=datastruct(1).monkey;
elseif length(varargin)==3
    datastruct=varargin{1};
    chamberWanted=varargin{2};
    nColumnsWanted=varargin{3};
    monkeyName=datastruct(1).monkey;
end
    % Discover all field names in the first element of datastruct
    %fields = fieldnames(datastruct(1));
    
    % Initialize chamberLeft and chamberRight with the same fields, but empty
    datastructLeft = [];%cell2struct(cell(length(fields), 1), fields);
    datastructRight = []; %cell2struct(cell(length(fields), 1), fields);
    
    leftIndex = 0; % Keep track of how many entries added to chamberLeft
    rightIndex = 0; % Keep track of how many entries added to chamberRight
    
switch monkeyName
    case {'Chip'}
        cutoffDate = '20230721';
        for entryNo = 1:length(datastruct)
            entry = datastruct(entryNo);
            entryDate = entry.date;
            if ~isempty(entryDate)
                if datenum(entryDate, 'yyyymmdd') >= datenum(cutoffDate, 'yyyymmdd')
                    % Check for chamberLeft criteria
                    if ~isempty(entry.powercycle) && ~isempty(entry.gammaCorrFactor) && ischar(entry.run) %~isempty(entry.gaussianContourLevel) && ~isempty(entry.powercycle) && ischar(entry.run)

                        % filtering for specific nColumns to analyze
                        if ~isempty(nColumnsWanted) && isfield(entry,'nColumnsWanted') && sum(entry.nColumnsWanted==nColumnsWanted)>=1
                            leftIndex = leftIndex + 1;
                            datastructLeft(leftIndex) = entryNo;
                        elseif isempty(nColumnsWanted)
                            leftIndex = leftIndex + 1;
                            datastructLeft(leftIndex) = entryNo;
                        end

                    end
                    excludeBlocks=[44 45 45 50 51 52:55 60 95];
                else
                    % Check for chamberRight criteria
                    if  ~isempty(entry.powercycle)
                        rightIndex = rightIndex + 1;
                        datastructRight(rightIndex) = entryNo;
                    end
                    excludeBlocks=10;
                end
            end
        end
    case {'Pepper'}
        cutoffDate = '20250404';
        for entryNo = 1:length(datastruct)
            entry = datastruct(entryNo);
            entryDate = entry.date;
            if ~isempty(entryDate)
                if datenum(entryDate, 'yyyymmdd') >= datenum(cutoffDate, 'yyyymmdd')
                    % Check for chamberLeft criteria
                    if ~isempty(entry.powercycle) && ~isempty(entry.gammaCorrFactor) && ischar(entry.run)  && strcmp(entry.chamber,chamberWanted) %~isempty(entry.gaussianContourLevel) && ~isempty(entry.powercycle) && ischar(entry.run)
                        % filtering for specific nColumns to analyze
                        if ~isempty(nColumnsWanted) && isfield(entry,'nColumnsWanted') && sum(entry.nColumnsWanted==nColumnsWanted)>=1
                            rightIndex = rightIndex + 1;
                            datastructRight(rightIndex) = entryNo;
                        elseif isempty(nColumnsWanted)
                            rightIndex = rightIndex + 1;
                            datastructRight(rightIndex) = entryNo;
                        end
                    end
                    excludeBlocks=nan;
                end
            end
        end
end

% Trim the unused preallocated entries
datastructLeft = datastructLeft(1:leftIndex);
datastructRight = datastructRight(1:rightIndex);

switch chamberWanted
    case 'L'
        datastructChamber=datastructLeft;
        analysisBlockID=setdiff(datastructChamber,excludeBlocks);
    case 'R'
        datastructChamber=datastructRight;
        analysisBlockID=setdiff(datastructChamber,excludeBlocks);
end
end
