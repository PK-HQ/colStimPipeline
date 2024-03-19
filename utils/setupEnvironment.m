function setupEnvironment()
%SETUPENVIRONMENT Configures the environment based on the operating system and computer name.

% Determine the operating system
if ispc
    % Get the computer name for Windows
    pcID = getenv('COMPUTERNAME');
    
    % Initialize userStr based on the computer ID
    switch pcID
        case {'CVIS-A64882','PSYC-A77304'}
            userStr = 'pktan';
        case {'LA-CPSD077020WD','LA-CPSA07019WD'}
            userStr = 'esexpt';
        otherwise
            error('Unknown PC ID: %s', pcID);
    end
    
    % Set directories for Windows
    monkeyDir = 'Y:/Chip/PsychometricData/';
    boxDataDir = ['C:/Users/' userStr '/Box/_Eyal/Columnar read-write/code/plots/Training/'];
    
    % Change to the monkeyDir directory
    cd(monkeyDir);
    
elseif ismac
    % Set directories for macOS
    monkeyDir = '~/Box/_Eyal/Columnar read-write/data/'; % Because data is here, not lab drive
    boxDataDir = '~/Box/_Eyal/Columnar read-write/code/plots/Training/';
    
    % Change to the monkeyDir directory
    cd(monkeyDir);
else
    error('Unsupported operating system.');
end

end
