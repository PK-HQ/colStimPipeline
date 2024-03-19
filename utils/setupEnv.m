function [mainPath,datastruct]=setupEnv(dataStructPath)
    % Determine the main path based on the system or hostname
    if ispc
        mainPath = 'Y:/';
    elseif contains(getenv('HOSTNAME'), 'psy.utexas.edu')
        mainPath = '/eslab/data/';
    else
        error('Unsupported system or hostname.');
    end

    % Construct the script path to be run
    scriptPath = fullfile(mainPath, dataStructPath);

    % Check if the script exists before attempting to run it
    if exist(scriptPath, 'file')
        run(scriptPath);
    else
        error('The specified script does not exist: %s', scriptPath);
    end

    % Misc
    set(0,'DefaultFigureWindowStyle','docked')
end
