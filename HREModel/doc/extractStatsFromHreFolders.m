function Tall = extractStatsFromHreFolders(rootDir, prefix, mode)
%extractStatsFromHreFolders Read [STATISTICS] logs from multiple HRE-bit folders.
%
% This function scans subfolders under rootDir, where each subfolder name is
% expected to be like:
%   0bit, 5bit, 10bit, ..., 70bit
% and each subfolder contains one or more *.log files for different SNR points.
%
% Usage:
%   Tall = extractStatsFromHreFolders(rootDir)
%   Tall = extractStatsFromHreFolders(rootDir, prefix)
%   Tall = extractStatsFromHreFolders(rootDir, prefix, mode)
%
% Inputs:
%   rootDir : parent folder that contains HRE subfolders (e.g., ./perf)
%   prefix  : passed to extractStatsFromLogs(); if empty, matches '*.log'
%   mode    : passed to extractStatsFromLogs(); 1=extended, 2=minimal
%
% Output:
%   Tall : table concatenating all rows across all HRE subfolders, with added
%          columns:
%            - HREbit (numeric)
%            - HREFolder (string)
%
% Notes / assumptions:
%   - HREbit is parsed from folder name using regexp '^(\d+)\s*bit$' (ignorecase).
%   - Parsing within each folder is delegated to doc/extractStatsFromLogs.m.

    if nargin < 1 || isempty(rootDir)
        rootDir = '.';
    end
    if nargin < 2
        prefix = '';
    end
    if nargin < 3 || isempty(mode)
        mode = 1;
    end

    if isstring(rootDir), rootDir = char(rootDir); end
    if isstring(prefix), prefix = char(prefix); end
    if isstring(mode), mode = str2double(char(mode)); end

    d = dir(rootDir);
    d = d([d.isdir]);
    d = d(~ismember({d.name}, {'.','..'}));

    allTables = {};
    for i = 1:numel(d)
        folderName = d(i).name;
        tok = regexp(folderName, '^(\d+)\s*bit$', 'tokens', 'once', 'ignorecase');
        if isempty(tok)
            continue;
        end

        hreBit = str2double(tok{1});
        subDir = fullfile(d(i).folder, folderName);

        try
            Tsub = extractStatsFromLogs(subDir, prefix, mode);
        catch ME
            % If a folder has no matching logs, skip it; otherwise rethrow.
            if strcmp(ME.identifier, 'extractStatsFromLogs:NoFiles')
                continue;
            end
            rethrow(ME);
        end

        Tsub.HREbit = repmat(hreBit, height(Tsub), 1);
        Tsub.HREFolder = repmat(string(folderName), height(Tsub), 1);
        allTables{end+1} = Tsub; %#ok<AGROW>
    end

    if isempty(allTables)
        error('extractStatsFromHreFolders:NoHreFolders', ...
            'No HRE subfolders matched the pattern ''<number>bit'' under: %s', rootDir);
    end

    Tall = vertcat(allTables{:});

    % Keep a stable column order (metadata first).
    metaCols = intersect({'HREbit','HREFolder','SNR','File','Path'}, Tall.Properties.VariableNames, 'stable');
    otherCols = setdiff(Tall.Properties.VariableNames, metaCols, 'stable');
    Tall = Tall(:, [metaCols, otherCols]);

    % Prefer sorting by HREbit then SNR then filename.
    if any(strcmp(Tall.Properties.VariableNames, 'SNR')) && any(~isnan(Tall.SNR))
        Tall = sortrows(Tall, {'HREbit','SNR','File'});
    else
        Tall = sortrows(Tall, {'HREbit','File'});
    end
end

