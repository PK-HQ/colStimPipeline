function displayTableWithHeaders(matrix, headers)
    % Check the size of the matrix
    if size(matrix, 2) ~= length(headers)
        error('Number of columns in the matrix must match the number of headers');
    end
    
    % Extract the first 5 rows of the first slice of the matrix
    data = matrix(1:5, :, 1);

    % Determine the maximum width for each column for proper alignment
    colWidths = max([cellfun(@length, headers); cellfun(@(x) length(num2str(x, '%.4f')), num2cell(data))], [], 1);
    
    % Create a format string for the header
    headerFormat = '|';
    for i = 1:length(headers)
        headerFormat = [headerFormat, ' %', num2str(colWidths(i)), 's |'];
    end
    headerFormat = [headerFormat, '\n'];
    
    % Create a format string for the data rows
    dataFormat = '|';
    for i = 1:size(data, 2)
        dataFormat = [dataFormat, ' %', num2str(colWidths(i)), '.4f |'];
    end
    dataFormat = [dataFormat, '\n'];
    
    % Print the header
    fprintf(headerFormat, headers{:});
    
    % Print a separator line
    separator = '|';
    for i = 1:length(colWidths)
        separator = [separator, repmat('-', 1, colWidths(i) + 2), '|'];
    end
    fprintf('%s\n', separator);
    
    % Print the data rows
    for i = 1:size(data, 1)
        fprintf(dataFormat, data(i, :));
    end
end
