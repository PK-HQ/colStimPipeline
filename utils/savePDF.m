function savePDF(filename, monkeyName, saveFlag, plotNo, plotTotal)

if saveFlag
    %% Saves plot as PDF
    % filename is a string
    % monkeyname is a string for monkey name
    % saveFlag == 1 means save
    % plotNo == 1 means save as a new file, >1 means append\
    
    if ispc
      mainPath='Y:/';
    elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
      mainPath='/eslab/data/';
    end
    
    %filenames
    filenameMonkeyPlot=[mainPath monkeyName '/Meta/' filename '-' num2str(plotNo) '.pdf'];
    filenameMonkey=[mainPath monkeyName '/Meta/' filename '.pdf'];
    filenameMeetings=[mainPath 'users/PK/Eyal/meetings/' filename '.pdf'];
    
    % Check dir exists, else make one
    checkDirExists(filenameMonkey)
    checkDirExists(filenameMeetings)

    %% Save as a new pdf per plot (if its not the final plot)
    export_fig(filenameMonkeyPlot,'-pdf','-nocrop'); offwarning()
    %fig=gcf; close(fig)

    %% If final plot, merge all plots that match filename and delete the individual pdfs
    if plotNo==plotTotal
        % Generate filenames
        for plotNo = 1:plotTotal+1
            filenameMonkeyPlots{plotNo} = [mainPath monkeyName '/Meta/' filename '-' num2str(plotNo) '.pdf'];
            filenameMonkeyPlots2{plotNo}=sprintf(filenameMonkeyPlots{plotNo},'\n');
        end
        % Merge into single pdf
        if exist(filenameMonkey,'file')
            delete(filenameMonkey)
        end
        append_pdfs(filenameMonkey,filenameMonkeyPlots{1:end})
        % copy to meetings folder
        copyfile(filenameMonkey, filenameMeetings)
        % Delete individual files
        if exist(filenameMonkey,'file') && exist(filenameMeetings,'file')
            for i = 1:length(filenameMonkeyPlots)-1
                delete(filenameMonkeyPlots{i}); % Delete each file
            end
        end
    end
end
end

function checkDirExists(fullPath)
    % This function extracts the directory path from a full file path, checks if 
    % this directory exists, and creates it if it does not.

    % Extract the directory path from the full path
    [folderPath, ~, ~] = fileparts(fullPath);

    % Check if the folder exists
    if exist(folderPath, 'dir') ~= 7
        % Folder does not exist, so create it
        mkdir(folderPath);
        fprintf('Directory created: %s\n', folderPath);
    else
        % Folder already exists
        %fprintf('Directory already exists: %s\n', folderPath);
    end
end
