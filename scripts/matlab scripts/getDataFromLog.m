%% Script to Process FER Simulation Log Files and Export to Excel
%
% This script reads simulation data from a structured directory, parses
% key statistics from .log files, and aggregates the results into a 
% single Excel spreadsheet.

% Author: Gemini
% Date:   2025-10-01

%% --- 1. Configuration ---
% Modify these variables to match your environment

% Set the base directory containing the 'matrixX' folders
baseDir = 'perf/perf_20x149'; 

% Set the desired name for the output Excel file
outputFileName = 'simulation_results.xlsx';

%% --- 2. Initialization ---
clc; % Clear command window
fprintf('Starting FER log processing...\n');

% Define the exact strings to search for within the log files
statKeys = {
    'Total simulated number', ...
    'RAW  BER', ...
    'LDPC BER', ...
    'LDPC FER', ...
    'MCRC FER', ...
    'DATA FER', ...
    'Retry Decoder average iterations'
};

% Initialize a cell array to store all the results
allResults = {};

%% --- 3. Find and Loop Through Matrix Folders ---

% Get a list of all subdirectories matching the 'matrix*' pattern
matrixFolders = dir(fullfile(baseDir, 'matrix*'));
matrixFolders = matrixFolders([matrixFolders.isdir]); % Keep only directories

if isempty(matrixFolders)
    warning('No directories matching "matrix*" found in "%s".', baseDir);
    return;
end

fprintf('Found %d matrix folders to process.\n', length(matrixFolders));

for i = 1:length(matrixFolders)
    currentMatrixFolder = matrixFolders(i).name;
    fprintf('\nProcessing: %s\n', currentMatrixFolder);
    
    % --- Extract Matrix ID from folder name using regular expressions
    matrixID_cell = regexp(currentMatrixFolder, '\d+', 'match');
    if isempty(matrixID_cell)
        fprintf('  -> Skipping folder "%s" (could not determine Matrix ID).\n', currentMatrixFolder);
        continue;
    end
    matrixID = str2double(matrixID_cell{1});
    
    % --- Find all .log files in the current matrix folder
    matrixPath = fullfile(baseDir, currentMatrixFolder);
    logFiles = dir(fullfile(matrixPath, 'snr_*.log'));
    
    if isempty(logFiles)
        fprintf('  -> No .log files found. Skipping this folder.\n');
        continue; % Skip to the next matrix folder
    end
    
    fprintf('  -> Found %d log files.\n', length(logFiles));
    
    %% --- 4. Loop Through Log Files and Parse Data ---
    for j = 1:length(logFiles)
        currentLogFile = logFiles(j).name;
        
        % --- Extract SNR value from file name
        snr_cell = regexp(currentLogFile, 'snr_(-?\d+\.?\d*)\.log', 'tokens');
        if isempty(snr_cell)
            fprintf('  -> Skipping file "%s" (could not parse SNR value).\n', currentLogFile);
            continue;
        end
        snrValue = str2double(snr_cell{1}{1});
        
        % --- Read and parse the contents of the log file
        logFilePath = fullfile(matrixPath, currentLogFile);
        try
            % Read all lines from the file (modern MATLAB function)
            fileLines = readlines(logFilePath); 
        catch
            fprintf('  -> ERROR: Could not read file "%s". Skipping.\n', logFilePath);
            continue;
        end
        
        % Find lines containing our statistics
        statsLines = fileLines(startsWith(fileLines, '[STATISTICS]'));
        
        % Prepare a temporary array to hold parsed values, default to NaN
        parsedValues = NaN(1, length(statKeys));
        
        for k = 1:length(statKeys)
            key = statKeys{k};
            % Find the line for the current key
            foundLine = statsLines(contains(statsLines, key));
            
            if ~isempty(foundLine)
                % Extract the numeric value after the colon
                value_cell = regexp(foundLine{1}, ':?\s*([0-9.eE+-]+)', 'tokens');
                if ~isempty(value_cell)
                    parsedValues(k) = str2double(value_cell{1}{1});
                end
            end
        end
        
        % --- Add the collected data to our main results cell array
        newRow = [{matrixID, snrValue}, num2cell(parsedValues)];
        allResults = [allResults; newRow];
    end
end

%% --- 5. Finalize and Write to Excel ---

if isempty(allResults)
    fprintf('\nProcessing complete. No data was found to write to Excel.\n');
    return;
end

% Convert the cell array to a MATLAB table for easy writing
try
    columnHeaders = {'Matrix_ID', 'SNR', 'Total_Simulated_Number', 'RAW_BER', ...
                     'LDPC_BER', 'LDPC_FER', 'MCRC_FER', 'DATA_FER', 'Avg_Iterations'};
    resultsTable = cell2table(allResults, 'VariableNames', columnHeaders);
    
    % Sort the table for better readability
    resultsTable = sortrows(resultsTable, {'Matrix_ID', 'SNR'});
    
    % Write the table to an Excel file
    writetable(resultsTable, outputFileName, 'Sheet', 1);
    
    fprintf('\n\n------------------------------------------------------\n');
    fprintf('Success! Data processing complete.\n');
    fprintf('Results exported to: %s\n', fullfile(pwd, outputFileName));
    fprintf('------------------------------------------------------\n');
    
catch ME
    fprintf('\n\nERROR: Could not write data to Excel file.\n');
    fprintf('Error message: %s\n', ME.message);
end