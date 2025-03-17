function formattedDates = getDates(daysBefore)
    % Ensure daysBefore is negative or zero
    if daysBefore > 0
        error('daysBefore should be a negative number or zero.');
    end

    % Get today's date
    today = datetime('today');

    % Calculate the start date based on the input daysBefore
    startDate = today + daysBefore;

    % Calculate the number of days in the range
    numDays = abs(daysBefore) + 1;

    % Initialize an array to store the formatted dates
    formattedDates = strings(numDays, 1);

    % Loop through each day from the start date to today
    for i = 0:abs(daysBefore)
        currentDate = startDate + days(i);
        formattedDates(i+1) = datestr(currentDate, 'yyyymmdd');
    end
end
