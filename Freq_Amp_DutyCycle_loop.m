% DUTY CYCLE, FREQUENCY, and AMPLITUDE
% This script calculates duty cycle as a proportion of time the cell is
% considered "ON" relative to the *entire* timecourse

% 07 04 2023: ADDED Amplitude and Frequency Analysis
% 02 19 2024: ADDED For Loops to Save Data

close all; clc;

% Load Timecourse
% directory where CalciumWaveForm is located

filepath = ('/Users/ibarevel/Desktop/Biocompatibility/');
filepath = ('//Users/ibarevel/Documents/MATLAB//Biocompatibility/Min6FunctionalityTest');
filepath = ('/Users/ibarevel/Documents/MATLAB//Biocompatibility/OutputMin6CaImaging');
javaaddpath('/Users/ibarevel/Documents/MATLAB/bioformats_package.jar');

%javaaddpath('/Users/ibarevel/Documents/MATLAB/estim/3-11-25/3-11-25Output/')
%filepath = ('D:/MyPassport/Evelyn_PhD/Estim/3-20-25/');
%load = ('/Users/ibarevel/Documents/MATLAB/estim/estim_output/CaWaveForm.mat');
%load('/Users/ibarevel/Documents/MATLAB/Biocompatibility/OutputDataMin6CaImaging/CaWaveForm.mat')

%load('/Users/ibarevel/Documents/MATLAB/Biocompatibility/OutputDataMin6CaImaging/CaWaveFormMed610_C_9-27-24.mat')

load('CaWaveForm.mat', 'CellTC')
% Assuming CellTC is your numeric array
CellTC(:, any(isnan(CellTC), 1)) = [];

%% Moving Average Smoothing
neighbors = 3; % Number of neighbors on each side
windowSize = 2 * neighbors + 1; % Total window size

% Apply moving average to each column
SmoothedCellTC = zeros(size(CellTC)); % Preallocate array
for i = 1:size(CellTC, 2)
    SmoothedCellTC(:, i) = movmean(CellTC(:, i), windowSize);
end

% Update CellTC with smoothed data
CellTC = SmoothedCellTC;

%% Define Dynamic Timeframe List
% You can add any list of timeframes for each run
timeframeList = { '1-130','131-310','311-1544', '1545-2161'}
%% Initialize Results Storage
numCells = size(CellTC, 2);
numTimeframes = length(timeframeList);

% Preallocate matrices to store results across timeframes
freqMatrix = NaN(numCells, numTimeframes);
ampMatrix = NaN(numCells, numTimeframes);
dutyCycle_M1_Matrix = NaN(numCells, numTimeframes);
dutyCycle_M2_Matrix = NaN(numCells, numTimeframes);

%% Loop through each timeframe in the list
for tf_idx = 1:numTimeframes
    % Extract start and end indices
    timeframe = timeframeList{tf_idx};
    timeframeParts = str2double(strsplit(timeframe, '-'));
    startIdx = timeframeParts(1);
    endIdx = timeframeParts(2);

    % Ensure valid indices
    if startIdx < 1 || endIdx > size(CellTC, 1)
        error('Timeframe indices are out of bounds.');
    end

    % Trim data for this timeframe
    CellTC_trimmed = CellTC(startIdx:endIdx, :);

    %% Duty Cycle Analysis (Method 1 & 2)
    for cell = 1:numCells
        Mean = CellTC_trimmed(:, cell);
        det = detrend(Mean);
        norm = normalize(det);
        avg_amp = mean(norm);
        on = avg_amp;

        % Find what percent of time the cell is "on"
        dutyCycle_M1_Matrix(cell, tf_idx) = sum(norm >= on) / length(Mean);

        % Method 2 using `dutycycle` function
        [d, ~, ~, ~, ~] = dutycycle(Mean);
        dutyCycle_M2_Matrix(cell, tf_idx) = mean(d, 'omitnan');
    end

    %% Frequency and Amplitude Analysis
    for cell = 1:numCells
        Mean = CellTC_trimmed(:, cell);
        s = detrend(Mean - mean(Mean));
        t = (0:length(s)-1);
        Fs = 1 / mean(diff(t));
        L = length(s);
        
        fts = fft(s);
        f = Fs * (0:(ceil(L/2))) / L;
        P2 = abs(fts / L);
        P1 = P2(1:(ceil(L/2) + 1));
        P1(2:end-1) = 2 * P1(2:end-1);

        highestFreq = f(find(P1 == max(P1), 1));
        value = find(f == highestFreq, 1);
        
        freqMatrix(cell, tf_idx) = f(value);
        ampMatrix(cell, tf_idx) = abs(fts(value));
    end
end






% % Step 1: Compute the correlation matrix
% corr_matrix = corr(CellTC);  % Compute pairwise correlation between cells
% 
% % Step 2: Plot the heatmap using the 'spring' colormap
% figure;
% heatmap(corr_matrix, 'Colormap', spring, 'ColorLimits', [-1 1]);
% title('Correlation Matrix of Cells: Entire Timeframe');
% xlabel('Cell');
% ylabel('Cell');
% 
% % Step 3: Save the figure as a JPG
% saveas(gcf, 'correlation_matrix_spring.jpg');  % Save the current figure as JPG
% 
% Assume CellTC is already defined
% 
% Assume CellTC is already defined with dimensions [num_timeframes x num_cells]
% CellTC: Rows are timeframes, columns are cells
% Define Dynamic Timeframe List
% Step 1: Compute the correlation matrix
corr_matrix = corr(CellTC);  % Compute pairwise correlation between cells

% Step 2: Compute the mean correlation of all cells
mean_corr = mean(corr_matrix(:));  % Flatten the correlation matrix and compute the mean

% Step 3: Plot the heatmap using the 'spring' colormap
figure;
h = heatmap(corr_matrix, 'Colormap', spring, 'ColorLimits', [-1 1]);
title('Correlation Matrix of Cells: Entire Timeframe');
xlabel('Cell');
ylabel('Cell');

% Step 4: Display the mean correlation as a separate annotation outside the heatmap
% Create an annotation box for the mean correlation
annotation('textbox', [0.5, 0.9, 0.1, 0.05], 'String', sprintf('Mean Correlation: %.2f', mean_corr), ...
           'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
           'FontWeight', 'bold', 'Color', 'black', 'EdgeColor', 'none');

% Step 5: Save the figure as a JPG
fig = figure('Renderer', 'painters');  % Use 2D software-based rendering

%saveas(gcf, 'correlation_matrix_spring.jpg');  % Save the current figure as JPG


% % Loop through each timeframe in the list
% for tf_idx = 1:numTimeframes
%     % Extract start and end indices for this timeframe
%     timeframe = timeframeList{tf_idx};
%     timeframeParts = str2double(strsplit(timeframe, '-'));
%     startIdx = timeframeParts(1);
%     endIdx = timeframeParts(2);
% 
%     % Ensure valid indices
%     if startIdx < 1 || endIdx > size(CellTC, 1)
%         error('Timeframe indices are out of bounds.');
%     end
% 
%     % Trim data for this timeframe (this will extract the rows for this timeframe)
%     CellTC_trimmed = CellTC(startIdx:endIdx, :);
% 
%     % Step 1: Compute the correlation matrix for this timeframe
%     corr_matrix = corr(CellTC_trimmed);  % Calculate correlation between cells for the timeframe
% 
%     % Step 2: Compute the mean correlation (excluding the diagonal)
%     % Extract the upper triangular part of the correlation matrix
%     triu_corr_matrix = triu(corr_matrix, 1);  % Exclude the diagonal (set to 0)
%     mean_corr = mean(triu_corr_matrix(triu_corr_matrix ~= 0));  % Calculate mean excluding zeros
% 
%     % Step 3: Plot the heatmap using the 'spring' colormap
%     figure;
%     heatmap(corr_matrix, 'Colormap', spring, 'ColorLimits', [-1 1]);
% 
%     % Add mean correlation value to the title
%     title(['Correlation Matrix of Cells (Timeframe ', num2str(tf_idx), ') - Mean Corr: ', num2str(mean_corr, '%.3f')]);
%     xlabel('Cell');
%     ylabel('Cell');
% 
%     % Step 4: Save the figure as a JPG for this timeframe
%     saveas(gcf, ['correlation_matrix_timeframe_', num2str(tf_idx), '.jpg']);
% 
%     % Close the figure for the next iteration (optional)
%    % close(gcf);
% end
% 
% 
% 
% %% Save DataOut File
% saveFile = input('Do you want to save the file? (Y/N): ', 's');
% 
% if strcmpi(saveFile, 'Y')
%     titleName = input('Enter the title for the file: ', 's');
%     fileName = sprintf('%s_DataOut.xlsx', titleName);
%     pathName = uigetdir('', 'Select Folder to Save the File');
% 
%     if isequal(pathName, 0)
%         disp('User canceled the file save.');
%     else
%         fullFileName = fullfile(pathName, fileName);
% 
%         try
%             % Convert matrices into tables with proper column labels
%             timeframeLabels = timeframeList;
% 
%             freqTable = array2table(freqMatrix, 'VariableNames', timeframeLabels);
%             ampTable = array2table(ampMatrix, 'VariableNames', timeframeLabels);
%             dutyCycleM1Table = array2table(dutyCycle_M1_Matrix, 'VariableNames', timeframeLabels);
%             dutyCycleM2Table = array2table(dutyCycle_M2_Matrix, 'VariableNames', timeframeLabels);
%         try
%             % Write tables to different sheets
%             writetable(freqTable, fullFileName, 'Sheet', 'Frequency');
%             writetable(ampTable, fullFileName, 'Sheet', 'Amplitude');
%             writetable(dutyCycleM1Table, fullFileName, 'Sheet', 'DutyCycleM1');
%             writetable(dutyCycleM2Table, fullFileName, 'Sheet', 'DutyCycleM2');
% 
%             % Check if CellTC is numeric or cell array
%             if isnumeric(CellTC)
%                 writematrix(CellTC, fullFileName, 'Sheet', 'CellTC'); % Use writematrix for numeric data
%             elseif iscell(CellTC)
%                 writecell(CellTC, fullFileName, 'Sheet', 'CellTC'); % Use writecell for cell arrays
%             else
%                 warning('CellTC is not a numeric matrix or cell array. Skipping save.');
%             end
% 
%             disp(['Analysis results saved to: ', fullFileName]);
%         catch ME
%             disp('Error while saving data to Excel:');
%             disp(ME.message);
        % end

%             % Write each result matrix to a different sheet
%             writetable(freqTable, fullFileName, 'Sheet', 'Frequency');
%             writetable(ampTable, fullFileName, 'Sheet', 'Amplitude');
%             writetable(dutyCycleM1Table, fullFileName, 'Sheet', 'DutyCycleM1');
%             writetable(dutyCycleM2Table, fullFileName, 'Sheet', 'DutyCycleM2');
%             writetable(CellTC, fullFileName, 'Sheet', 'CellTC');
%             disp(['Analysis results saved to: ', fullFileName]);
%         catch ME
%             disp('Error while saving data to Excel:');
%             disp(ME.message);
     %   end
     % end
% elseif strcmpi(saveFile, 'N')
%     disp('File save operation canceled. Exiting.');
% else
%     disp('Invalid input. Exiting.');
 % end
% % Loop through each timeframe in the list
% meanCorrValues = zeros(numTimeframes, 1); % Store mean correlation values
% timeframeLabels = cell(numTimeframes, 1); % Store timeframe labels

for tf_idx = 1:numTimeframes
    % Extract start and end indices for this timeframe
    timeframe = timeframeList{tf_idx};
    timeframeParts = str2double(strsplit(timeframe, '-'));
    startIdx = timeframeParts(1);
    endIdx = timeframeParts(2);

    % Ensure valid indices
    if startIdx < 1 || endIdx > size(CellTC, 1)
        error('Timeframe indices are out of bounds.');
    end

    % Trim data for this timeframe
    CellTC_trimmed = CellTC(startIdx:endIdx, :);

    % Step 1: Compute the correlation matrix
    corr_matrix = corr(CellTC_trimmed);

    % Step 2: Compute the mean correlation (excluding diagonal)
    triu_corr_matrix = triu(corr_matrix, 1);
    mean_corr = mean(triu_corr_matrix(triu_corr_matrix ~= 0));

    % Store mean correlation value
    meanCorrValues(tf_idx) = mean_corr;
    timeframeLabels{tf_idx} = timeframe;

    % Step 3: Plot using imagesc (better for exportgraphics)
    fig = figure('Renderer', 'painters');
    imagesc(corr_matrix);
    colorbar;
    colormap(spring);
    caxis([-1 1]);
    axis equal tight;
    
    title(['Correlation Matrix (Timeframe ', num2str(tf_idx), ') - Mean Corr: ', num2str(mean_corr, '%.3f')]);
    xlabel('Cell');
    ylabel('Cell');

    ax = gca;
    disableDefaultInteractivity(ax);  % Prevent GUI artifacts
    set(ax, 'XTick', 1:size(corr_matrix,1), ...
            'YTick', 1:size(corr_matrix,2), ...
            'XTickLabelRotation', 90);

    % Step 4: Save figure
    exportgraphics(fig, ['correlation_matrix_timeframe_', num2str(tf_idx), '.jpg'], 'Resolution', 300);
    close(fig);  % Optional: close figure to save memory
end


%% Save DataOut File
saveFile = input('Do you want to save the file? (Y/N): ', 's');

if strcmpi(saveFile, 'Y')
    titleName = input('Enter the title for the file: ', 's');
    fileName = sprintf('%s_DataOut.xlsx', titleName);
    pathName = uigetdir('', 'Select Folder to Save the File');

    if isequal(pathName, 0)
        disp('User canceled the file save.');
    else
        fullFileName = fullfile(pathName, fileName);

        try
            % Convert matrices into tables
            % Before creating freqTable
            nCols = size(freqMatrix, 2);
            
            if ~exist('timeframeLabels', 'var')
                warning('Variable timeframeLabels does not exist. Generating default variable names.');
                autoLabels = strcat("T", string(1:nCols));
            elseif length(timeframeLabels) ~= nCols
                warning('Mismatch detected: Number of timeframeLabels (%d) does not match number of columns (%d).', length(timeframeLabels), nCols);
                disp('timeframeLabels content:');
                disp(timeframeLabels);
                disp('Auto-generating default variable names.');
                autoLabels = strcat("T", string(1:nCols));
            else
                autoLabels = string(timeframeLabels);
            end
            
            freqTable = array2table(freqMatrix, 'VariableNames', autoLabels);
            
            % Repeat similar checks for ampMatrix, dutyCycle_M1_Matrix, dutyCycle_M2_Matrix...

            ampTable = array2table(ampMatrix, 'VariableNames', timeframeLabels);
            dutyCycleM1Table = array2table(dutyCycle_M1_Matrix, 'VariableNames', timeframeLabels);
            dutyCycleM2Table = array2table(dutyCycle_M2_Matrix, 'VariableNames', timeframeLabels);
            
            % Write tables to different sheets
            writetable(freqTable, fullFileName, 'Sheet', 'Frequency');
            writetable(ampTable, fullFileName, 'Sheet', 'Amplitude');
            writetable(dutyCycleM1Table, fullFileName, 'Sheet', 'DutyCycleM1');
            writetable(dutyCycleM2Table, fullFileName, 'Sheet', 'DutyCycleM2');

            % Check if CellTC is numeric or cell array
            if isnumeric(CellTC)
                writematrix(CellTC, fullFileName, 'Sheet', 'CellTC'); 
            elseif iscell(CellTC)
                writecell(CellTC, fullFileName, 'Sheet', 'CellTC'); 
            else
                warning('CellTC is not a numeric matrix or cell array. Skipping save.');
            end
            
            % Save correlation matrices
            writematrix(meanCorrValues, fullFileName, 'Sheet', 'Mean_Correlation', 'Range', 'B2');
            writecell(timeframeLabels, fullFileName, 'Sheet', 'Mean_Correlation', 'Range', 'A2');
            writecell({'Timeframe', 'Mean Pearson Correlation'}, fullFileName, 'Sheet', 'Mean_Correlation', 'Range', 'A1');

            disp(['Analysis results saved to: ', fullFileName]);
        catch ME
            disp('Error while saving data to Excel:');
            disp(ME.message);
        end
    end
end
