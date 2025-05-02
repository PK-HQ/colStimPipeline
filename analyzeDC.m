%load('Y:\Pepper\Pepper20250417\run1\M32D20250417R1StabDFFTAmpS004E023PF0400.mat')
%load('Y:\Pepper\Pepper20250417\run1\M32D20250417R1TS.mat')

%% Define conds from TS file
blankID=find(TS.Header.Conditions.TypeCond==0);
condID=find(TS.Header.Conditions.TypeCond>0);
bmpData=extractParams(TS.Header.Conditions.ProjImg);
phs0=find(bmpData.Phs(condID)==0);
phs90=find(bmpData.Phs(condID)==180);
DCs=bmpData.DCs(condID);
DCs=DCs(1:numel(DCs)/2);
%% Process imaging data
% Subtract blank
condData=DataCond(:,:,condID)-DataCond(:,:,blankID);
% Smooth the map using nanconv (handles NaNs)
filtWidth = filtSigma*6;
imageFilter = fspecial('gaussian',filtWidth,filtSigma);
for i=1:numel(condID)
    condData(:,:,i) = nanconv(condData(:,:,i), imageFilter, 'nanout'); % USING NANCONV
end


% F0 (sum of phases)
F0=condData(:,:,phs0)+condData(:,:,phs90);
% F1 (difference of phases)
F1=condData(:,:,phs90)-condData(:,:,phs0);
% F1/F0
F1F0=F1./F0;
F1F0mean=squeeze(mean(F1F0,[1 2]));



%% Plotting Functionality
figure; % Create a new figure window
set(gcf, 'Name', 'F0, F1, F1/F0 Maps', 'NumberTitle', 'off'); % Set figure title

data_to_plot = {F0, F1, F1F0}; % Cell array holding the data matrices
row_labels = {'F0', 'F1', 'F1/F0'}; % Labels for rows/y-axis
num_rows = 3;

% Determine number of columns based on the data
if size(F0, 3) == numel(DCs) && size(F1, 3) == numel(DCs) && size(F1F0, 3) == numel(DCs)
    num_cols = numel(DCs);
    if num_cols ~= 5
       warning('Number of DCs (%d) does not match requested 5 columns. Adjusting grid.', num_cols);
       % You might want to add logic here to handle cases where num_cols is not 5,
       % like plotting only the first 5, or adjusting the subplot grid.
       % For now, we'll proceed assuming num_cols is the intended number.
    end
else
    error('The third dimension of F0, F1, F1F0 must match the number of unique DCs.');
end


for r = 1:num_rows % Loop through rows (F0, F1, F1F0)
    current_data = data_to_plot{r};

    % Determine reasonable, consistent color limits for the current row
    finite_data = current_data(isfinite(current_data(:))); % Exclude NaN/Inf
    if isempty(finite_data)
         min_val = 0;
         max_val = 1;
    else
        min_val = min(finite_data);
        max_val = max(finite_data);
    end

    if r == 3 % F1/F0 often needs symmetric limits around 0
        max_abs = max(abs([min_val, max_val]));
        if max_abs == 0 % Handle case where all values are zero
             c_limits = [-1, 1]; % Assign default limits
        else
             c_limits = [-max_abs, max_abs];
        end
    else % F0, F1 limits (adjust if negative values are expected)
        % Defaulting to full range; adjust if needed e.g., c_limits = [0, max_val]
        if min_val == max_val % Handle case where all values are the same
            c_limits = [min_val - 0.5, max_val + 0.5]; % Add padding
             if min_val == 0
                 c_limits = [-0.5, 0.5];
             end
        else
             c_limits = [min_val, max_val];
        end
    end
     % Ensure limits are not NaN and min < max
     if any(isnan(c_limits)) || c_limits(1) >= c_limits(2)
         c_limits = [0 1]; % Fallback default
     end


    for c = 1:num_cols % Loop through columns (DCs)
        plot_index = (r - 1) * num_cols + c;
        subplot(num_rows, num_cols, plot_index);

        imagesc(current_data(:, :, c));
        axis image; % Correct aspect ratio for images
        colorbar;
        try % Set color limits, handle potential errors if limits are invalid
           caxis(c_limits);
        catch ME
           warning('Could not set caxis for subplot (%d, %d): %s', r, c, ME.message);
        end


        % Add titles only for the first row (r=1)
        if r == 1
            title(sprintf('DC = %.2f', DCs(c)));
        end

        % Add y-axis labels only for the first column (c=1)
        if c == 1
            ylabel(row_labels{r});
        end

        % Optional: Remove axis ticks for a cleaner image grid
        set(gca, 'XTick', [], 'YTick', []);

        % Optional: Set a colormap (e.g., for F1/F0)
        if r == 3
            colormap(gca, 'cool'); % Example: using 'cool' colormap for F1/F0
        else
             colormap(gca, 'parula'); % Default or another map for F0/F1
        end
    end
end

% Optional: Add an overall title to the figure
% Adjust position: [left, bottom, width, height] in normalized units
% sgtitle('F0, F1, and F1/F0 Maps across Duty Cycles');
% Or use annotation for better control:
annotation('textbox', [0 0.9 1 0.1], ...
           'String', 'F0, F1, and F1/F0 Maps across Duty Cycles', ...
           'EdgeColor', 'none', ...
           'HorizontalAlignment', 'center', ...
           'FontSize', 14, ...
           'FontWeight', 'bold');







%% Subfuncs
function bmpData=extractParams(proj_img_list)
num_files = numel(proj_img_list);

% Initialize arrays to store extracted data
% Using NaN as a placeholder if a pattern is not found
bmpData.DCs = NaN(num_files, 1);
bmpData.SFs = NaN(num_files, 1);
bmpData.Phs = NaN(num_files, 1);
bmpData.referenceDCs = NaN(num_files, 1);
bmpData.bmpFilenames = proj_img_list(:); % Store filenames directly (as column cell array)

% Define regular expression patterns
sf_pattern = 'CBS(\d{3,})';     % Capture 3 or more digits after 'CBS'
phs_pattern = 'P(\d{5,})';     % Capture 3 or more digits after 'CBS'
dc_pattern = 'D(\d{3,})L';       % Capture 3 or more digits after 'D' and before 'L'
ref_dc_pattern = '_DC(\d+)ref';  % Capture digits after '_DC' and before 'ref'

% Iterate through each filename string in the cell array
for i = 1:num_files
    filename = proj_img_list{i};

    % Extract SF
    tokens = regexp(filename, sf_pattern, 'tokens', 'once');
    if ~isempty(tokens)
        bmpData.SFs(i) = str2double(tokens{1}) / 100.0;
    end
    
    % Extract Phase
    tokens = regexp(filename, phs_pattern, 'tokens', 'once');
    if ~isempty(tokens)
        bmpData.Phs(i) = str2double(tokens{1}) / 100.0;
    end

    % Extract DC
    tokens = regexp(filename, dc_pattern, 'tokens', 'once');
    if ~isempty(tokens)
        bmpData.DCs(i) = str2double(tokens{1}) / 100.0;
    end

    % Extract referenceDC
    tokens = regexp(filename, ref_dc_pattern, 'tokens', 'once');
    if ~isempty(tokens)
        bmpData.referenceDCs(i) = str2double(tokens{1}) / 100.0;
    end
end
end