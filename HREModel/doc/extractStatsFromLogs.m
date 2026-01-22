function T = extractStatsFromLogs(logDir, prefix, mode)
%extractStatsFromLogs Extract selected [STATISTICS] metrics from log files.
%
% Usage:
%   T = extractStatsFromLogs(logDir, prefix)
%   T = extractStatsFromLogs(logDir, prefix, mode)
%
% Inputs:
%   logDir : folder containing log files
%   prefix : file prefix, e.g. 'test' will match 'test*.log'
%            if empty, will match '*.log'
%   mode   : parsing mode
%            1 = current (extended) format (default)
%            2 = minimal format (only a small set of fields)
%
% Output:
%   T : table, one row per log file.
%
%   mode=1 includes:
%       RAW_BER, TheoRBER, FBC_total_x0_x2, FBC_TAWGN_x1_x2, HRE_like,
%       FixedHRE_flip_n, BINSTAT0 (0->7/5/3/1/self), BINSTAT1 (1->0/6/4/2/self),
%       LDPC_FER, RetryDecAvgIter.
%
%   mode=2 includes only:
%       RAW_BER, TheoRBER, FBC_bit, FBC_noHRE, HRE_like, FixedHRE, LDPC_FER.
%
% Notes:
%   - This parser searches the last occurrence of each [STATISTICS] field.
%   - If [STATISTICS] is missing (e.g., simulation did not finish), it falls
%     back to the last *complete* [SIM] summary block.
%   - It also tries to infer SNR from filename like '*snr3.7*.log'; if not
%     found, it tries to parse SNR from the log content; otherwise SNR=NaN.

    if nargin < 1 || isempty(logDir)
        logDir = '.';
    end
    if nargin < 2
        prefix = '';
    end
    if nargin < 3 || isempty(mode)
        mode = 1;
    end

    if isstring(logDir), logDir = char(logDir); end
    if isstring(prefix), prefix = char(prefix); end
    if isstring(mode), mode = str2double(char(mode)); end
    if ~isscalar(mode) || ~ismember(mode, [1 2])
        error('extractStatsFromLogs:BadMode', 'mode must be 1 or 2.');
    end

    if isempty(prefix)
        fileGlob = '*.log';
    else
        fileGlob = [prefix, '*.log'];
    end

    files = dir(fullfile(logDir, fileGlob));
    if isempty(files)
        error('extractStatsFromLogs:NoFiles', 'No log files matched: %s', fullfile(logDir, fileGlob));
    end

    rows = repmat(makeEmptyRow(mode), 0, 1);
    for i = 1:numel(files)
        filePath = fullfile(files(i).folder, files(i).name);
        txt = fileread(filePath);
        row = parseOne(txt, files(i).name, filePath, mode);
        rows(end+1,1) = row; %#ok<AGROW>
    end

    T = struct2table(rows);

    % Prefer sorting by SNR if available; otherwise by filename.
    if any(strcmp(T.Properties.VariableNames, 'SNR')) && any(~isnan(T.SNR))
        T = sortrows(T, {'SNR','File'});
    else
        T = sortrows(T, 'File');
    end
end

function row = makeEmptyRow(mode)
    common = struct( ...
        'File', "", ...
        'Path', "", ...
        'SNR', NaN, ...
        'RAW_BER', NaN, ...
        'TheoRBER', NaN, ...
        'HRE_like', NaN, ...
        'LDPC_FER', NaN, ...
        'RetryDecAvgIter', NaN);

    if mode == 2
        row = common;
        row.FBC_bit = NaN;
        row.FBC_noHRE = NaN;
        row.FixedHRE = NaN;
        return;
    end

    % mode == 1 (extended)
    row = common;
    row.FBC_total_x0_x2 = NaN;
    row.FBC_TAWGN_x1_x2 = NaN;
    row.FixedHRE_flip_n = NaN;
    row.BIN0_to7_HRE = NaN;
    row.BIN0_to5 = NaN;
    row.BIN0_to3 = NaN;
    row.BIN0_to1 = NaN;
    row.BIN0_self_0246 = NaN;
    row.BIN1_to0_HRE = NaN;
    row.BIN1_to6 = NaN;
    row.BIN1_to4 = NaN;
    row.BIN1_to2 = NaN;
    row.BIN1_self_1357 = NaN;
end

function row = parseOne(txt, fileName, filePath, mode)
    row = makeEmptyRow(mode);
    row.File = string(fileName);
    row.Path = string(filePath);
    row.SNR = parseSNR(fileName, txt);

    prefix = '^\s*\[STATISTICS\]\s+';

    row.RAW_BER = lastNumber(txt, [prefix 'RAW\s+BER\s*:\s*([0-9eE+\-\.]+)']);
    row.TheoRBER = lastNumber(txt, [prefix 'TheoRBER\s*:\s*([0-9eE+\-\.]+)']);
    row.HRE_like = lastNumber(txt, [prefix 'HRE-like\s*:\s*([0-9eE+\-\.]+)']);
    row.LDPC_FER = lastNumber(txt, [prefix 'LDPC\s+FER\s*:\s*([0-9eE+\-\.]+)']);
    row.RetryDecAvgIter = lastNumber(txt, [prefix 'Retry decoder average iteration\s*:\s*([0-9eE+\-\.]+)']);

    if mode == 2
        row.FBC_bit = lastNumber(txt, [prefix 'FBC\(bit\)\s*:\s*([0-9eE+\-\.]+)']);
        row.FBC_noHRE = lastNumber(txt, [prefix 'FBC\(no\s+HRE\)\s*:\s*([0-9eE+\-\.]+)']);
        row.FixedHRE = lastNumber(txt, [prefix 'Fixed\s+HRE\s*:\s*([0-9eE+\-\.]+)']);
        row = fillFromSimIfMissing(row, txt, mode);
        return;
    end

    % mode == 1 (extended)
    row.FBC_total_x0_x2 = lastNumber(txt, [prefix 'FBC\(total,x0->x2\)\s*:\s*([0-9eE+\-\.]+)']);
    row.FBC_TAWGN_x1_x2 = lastNumber(txt, [prefix 'FBC\(TAWGN,x1->x2\)\s*:\s*([0-9eE+\-\.]+)']);
    row.FixedHRE_flip_n = lastNumber(txt, [prefix 'Fixed\s+HRE\s*\(flip\s+n\)\s*:\s*([0-9eE+\-\.]+)']);

    % BINSTAT0
    tok0 = lastTokens(txt, [prefix 'BINSTAT0:\s*' ...
        '0->7\(HRE\)\s+([0-9eE+\-\.]+)\s*\|\s*' ...
        '0->5\s+([0-9eE+\-\.]+)\s*\|\s*' ...
        '0->3\s+([0-9eE+\-\.]+)\s*\|\s*' ...
        '0->1\s+([0-9eE+\-\.]+)\s*\|\s*' ...
        '0->self\(0/2/4/6\)\s+([0-9eE+\-\.]+)']);
    if ~isempty(tok0)
        row.BIN0_to7_HRE = str2double(tok0{1});
        row.BIN0_to5 = str2double(tok0{2});
        row.BIN0_to3 = str2double(tok0{3});
        row.BIN0_to1 = str2double(tok0{4});
        row.BIN0_self_0246 = str2double(tok0{5});
    end

    % BINSTAT1
    tok1 = lastTokens(txt, [prefix 'BINSTAT1:\s*' ...
        '1->0\(HRE\)\s+([0-9eE+\-\.]+)\s*\|\s*' ...
        '1->6\s+([0-9eE+\-\.]+)\s*\|\s*' ...
        '1->4\s+([0-9eE+\-\.]+)\s*\|\s*' ...
        '1->2\s+([0-9eE+\-\.]+)\s*\|\s*' ...
        '1->self\(1/3/5/7\)\s+([0-9eE+\-\.]+)']);
    if ~isempty(tok1)
        row.BIN1_to0_HRE = str2double(tok1{1});
        row.BIN1_to6 = str2double(tok1{2});
        row.BIN1_to4 = str2double(tok1{3});
        row.BIN1_to2 = str2double(tok1{4});
        row.BIN1_self_1357 = str2double(tok1{5});
    end

    row = fillFromSimIfMissing(row, txt, mode);
end

function v = lastNumber(txt, pattern)
% Return the last numeric token matched by pattern; NaN if not found.
    tokens = regexp(txt, pattern, 'tokens', 'lineanchors');
    if isempty(tokens)
        v = NaN;
        return;
    end
    v = str2double(tokens{end}{1});
end

function tok = lastTokens(txt, pattern)
% Return the last token list matched by pattern; {} if not found.
    tokens = regexp(txt, pattern, 'tokens', 'lineanchors');
    if isempty(tokens)
        tok = {};
        return;
    end
    tok = tokens{end};
end

function snr = parseSNR(fileName, txt)
% Try filename first (e.g. *_snr3.7.log), then CH_TRX line.
    snr = NaN;

    tok = regexp(fileName, 'snr([-+]?\d+(?:\.\d+)?)', 'tokens', 'once', 'ignorecase');
    if ~isempty(tok)
        snr = str2double(tok{1});
        return;
    end

    tok = regexp(txt, '\[CH_TRX\].*SNR\s*=\s*([-+]?\d+(?:\.\d+)?)', 'tokens', 'once');
    if ~isempty(tok)
        snr = str2double(tok{1});
        return;
    end

    tok = regexp(txt, '^\[SIM\]\s+SNR\s*:\s*([-+]?\d+(?:\.\d+)?)', 'tokens', 'once', 'lineanchors');
    if ~isempty(tok)
        snr = str2double(tok{1});
    end
end

function row = fillFromSimIfMissing(row, txt, mode)
% Fill missing fields from the last complete [SIM] block (if any).
    sim = lastCompleteSimBlock(txt, mode);
    if isempty(sim)
        return;
    end

    names = fieldnames(sim);
    for i = 1:numel(names)
        f = names{i};
        if isfield(row, f) && isnan(row.(f)) && ~isnan(sim.(f))
            row.(f) = sim.(f);
        end
    end
end

function sim = lastCompleteSimBlock(txt, mode)
% Return a struct of fields parsed from the last *complete* [SIM] summary block.
% If none found, return [].
%
% We define a SIM block by lines starting at:
%   [SIM] Statistical result of ...
% and ending right before the next such line or EOF.

    sim = [];
    starts = regexp(txt, '^\s*\[SIM\]\s+Statistical result of .*$', 'start', 'lineanchors');
    if isempty(starts)
        return;
    end
    starts = starts(:);

    ends = [starts(2:end)-1; numel(txt)];
    for k = numel(starts):-1:1
        blk = txt(starts(k):ends(k));
        cand = parseSimFields(blk, mode);
        if isSimComplete(cand, mode)
            sim = cand;
            return;
        end
    end
end

function s = parseSimFields(blk, mode)
% Parse required fields from a SIM block; missing fields become NaN.
    s = struct();
    sp = '^\s*\[SIM\]\s+';

    s.RAW_BER = lastNumber(blk, [sp 'RAW\s+BER\s*:\s*([0-9eE+\-\.]+)']);
    s.TheoRBER = lastNumber(blk, [sp 'TheoRBER\s*:\s*([0-9eE+\-\.]+)']);
    s.HRE_like = lastNumber(blk, [sp 'HRE-like\s*:\s*([0-9eE+\-\.]+)']);
    s.LDPC_FER = lastNumber(blk, [sp 'LDPC\s+FER\s*:\s*([0-9eE+\-\.]+)']);
    s.RetryDecAvgIter = lastNumber(blk, [sp 'Retry decoder average iteration\s*:\s*([0-9eE+\-\.]+)']);

    if mode == 2
        s.FBC_bit = lastNumber(blk, [sp 'FBC\(bit\)\s*:\s*([0-9eE+\-\.]+)']);
        s.FBC_noHRE = lastNumber(blk, [sp 'FBC\(no\s+HRE\)\s*:\s*([0-9eE+\-\.]+)']);
        s.FixedHRE = lastNumber(blk, [sp 'Fixed\s+HRE\s*:\s*([0-9eE+\-\.]+)']);
        return;
    end

    s.FBC_total_x0_x2 = lastNumber(blk, [sp 'FBC\(total,x0->x2\)\s*:\s*([0-9eE+\-\.]+)']);
    s.FBC_TAWGN_x1_x2 = lastNumber(blk, [sp 'FBC\(TAWGN,x1->x2\)\s*:\s*([0-9eE+\-\.]+)']);
    s.FixedHRE_flip_n = lastNumber(blk, [sp 'Fixed\s+HRE\s*\(flip\s+n\)\s*:\s*([0-9eE+\-\.]+)']);
end

function tf = isSimComplete(s, mode)
% Determine whether the SIM block contains a complete set of required fields.
    if mode == 2
        req = {'RAW_BER','TheoRBER','FBC_bit','FBC_noHRE','HRE_like','FixedHRE','LDPC_FER'};
    else
        req = {'RAW_BER','TheoRBER','FBC_total_x0_x2','FBC_TAWGN_x1_x2','HRE_like','FixedHRE_flip_n','LDPC_FER'};
    end
    tf = true;
    for i = 1:numel(req)
        f = req{i};
        if ~isfield(s, f) || isnan(s.(f))
            tf = false;
            return;
        end
    end
end
