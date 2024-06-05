function table2excel(data_table, filename)
    % Function to export a MATLAB table to an Excel file
    % Arguments:
    %   data_table: A MATLAB table that you want to export
    %   filename: Name of the Excel file to create
    
    % Check if the filename ends with '.xlsx', append if not
    if ~endsWith(filename, '.xlsx')
        filename = [filename '.xlsx'];
    end
    
    % Export the table to an Excel file
    try
        writetable(data_table, filename);
        fprintf('Table exported successfully to %s\n', filename);
    catch ME
        fprintf('Failed to export table: %s\n', ME.message);
    end
end
