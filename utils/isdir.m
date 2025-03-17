function result = isdir(folderPath)
    % ISDIR  Checks if the given path is a directory, made to allow
    % backwards compatibility with isfolder in new matlab versions
    %   RESULT = ISDIR(FOLDERPATH) returns true if FOLDERPATH is a valid directory,
    %   and false otherwise. This function replicates the behavior of the removed
    %   `isdir` function by using `isfolder`.

    % Check if the input path is a folder
    result = isfolder(folderPath);
end
