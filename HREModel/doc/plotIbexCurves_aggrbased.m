function plotIbexCurves_aggrbased(ferSummaryXlsx, tpSummaryXlsx, varargin)
%plotIbexCurves_aggrbased  Plot FER/AvgIter/Throughput curves from aggregate_results.py summary xlsx.
%
% Inputs:
%   ferSummaryXlsx: path to FER_summary_*.xlsx (or a directory containing it)
%   tpSummaryXlsx : path to Throughput_summary_*.xlsx (or a directory containing it; can be "" to skip TP plots)
%
% The summary xlsx format is produced by aggregate_results.py:
%   - One sheet per code rate (e.g., "1952B"), plus a "Summary" sheet.
%   - Each code-rate sheet contains repeated blocks:
%       "Code Rate: ... | Decoder: ... | Channel: ..." (header row)
%       metric rows: RAW_BER, LDPC_FER, FAIL_CW, AvgIter, ...
%     where each metric row contains all test points across columns.
%
% Example:
%   plotIbexCurves_aggrbased("path/summary/FER_summary_20250101_000000.xlsx", ...
%                            "path/summary/Throughput_summary_20250101_000000.xlsx", ...
%                            "MinFailCW", 4, ...
%                            "FigurePrefix", "IBEX_");

%% Parse options
p = inputParser;
p.addParameter("Model", "AUTO"); % "AUTO" | "AWGN" | "ERR_INJ"
p.addParameter("CodeRates", {}, @(x) iscell(x) || isstring(x));
p.addParameter("Configs", {}, @(x) iscell(x) || isstring(x));
p.addParameter("MinFailCW", 4, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter("KDirection", "increasing", @(s) any(strcmpi(string(s), ["increasing","decreasing"])));
p.addParameter("FigurePrefix", "", @(s) ischar(s) || isstring(s));
p.addParameter("Verbose", false, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;
figPrefix = string(opt.FigurePrefix);

ferSummaryXlsx = resolveSummaryXlsx(string(ferSummaryXlsx), "FER_summary_*.xlsx");
tpSummaryXlsx = resolveSummaryXlsx(string(tpSummaryXlsx), "Throughput_summary_*.xlsx");

%% Load FER blocks from summary workbook
[ratesAll, blocksByRate] = loadSummaryWorkbookBlocks(ferSummaryXlsx);
[rates, cfgs, ferData] = buildFerDataFromBlocks(ratesAll, blocksByRate, opt.CodeRates, opt.Configs, opt.MinFailCW, opt.Model);
if isempty(rates) || isempty(cfgs)
    error("plotIbexCurves_aggrbased:NoData", "No usable FER data found in %s.", ferSummaryXlsx);
end
if opt.Verbose
    fprintf("[plotIbexCurves_aggrbased] FER_summary=%s\n", ferSummaryXlsx);
    fprintf("[plotIbexCurves_aggrbased] rates=%d, configs=%d\n", numel(rates), numel(cfgs));
    fprintf("[plotIbexCurves_aggrbased] configs: %s\n", strjoin(cfgs, ", "));
end

%% Plot 1: RBER-FER (AWGN/TAWGN) or K-FER (ERR_INJ)
for ci = 1:numel(ferData)
    cfg = ferData(ci).config;
    if isempty(ferData(ci).rates)
        continue;
    end

    [model, axisName] = inferModelAndAxis(opt.Model, ferData(ci).rates(1).T);

    f = figure("Name", figPrefix + "FER_" + cfg, "NumberTitle", "off", "Color", "w");
    ax = axes(f);
    hold(ax, "on");
    grid(ax, "on");

    leg = strings(0, 1);
    for ri = 1:numel(ferData(ci).rates)
        rate = ferData(ci).rates(ri).rate;
        T = ferData(ci).rates(ri).T;

        [x, y] = selectFerXY(T, model);
        if isempty(x)
            continue;
        end
        [x, y] = sortByX(x, y);

        semilogy(ax, x, y, "-o", "LineWidth", 1.2, "MarkerSize", 4);
        leg(end+1) = rate; %#ok<AGROW>
    end

    xlabel(ax, axisName, "Interpreter", "none");
    ylabel(ax, "FER (LDPC_FER)", "Interpreter", "none");
    title(ax, figPrefix + cfg + " | " + model, "Interpreter", "none");
    set(ax, "YScale", "log");
    if strcmpi(model, "ERR_INJ") && strcmpi(opt.KDirection, "decreasing")
        set(ax, "XDir", "reverse");
    end
    if ~isempty(leg)
        legend(ax, leg, "Interpreter", "none", "Location", "best");
    end
end

%% Load throughput blocks (optional)
tpData = struct();
if strlength(tpSummaryXlsx) > 0 && isfile(tpSummaryXlsx)
    [tpRatesAll, tpBlocksByRate] = loadSummaryWorkbookBlocks(tpSummaryXlsx);
    [~, ~, tpData] = buildTpDataFromBlocks(tpRatesAll, tpBlocksByRate, opt.CodeRates, opt.Configs);
end
if opt.Verbose && strlength(tpSummaryXlsx) > 0 && isfile(tpSummaryXlsx)
    fprintf("[plotIbexCurves_aggrbased] Throughput_summary=%s\n", tpSummaryXlsx);
end

%% Plot 2: RBER/K - AvgIter (from Throughput_summary blocks)
for ci = 1:numel(tpData)
    cfg = tpData(ci).config;
    if isempty(tpData(ci).rates)
        continue;
    end

    [model, axisName] = inferModelFromTp(tpData(ci).rates(1).Ttp);

    f = figure("Name", figPrefix + "AvgIter_" + cfg, "NumberTitle", "off", "Color", "w");
    ax = axes(f);
    hold(ax, "on");
    grid(ax, "on");

    leg = strings(0, 1);
    for ri = 1:numel(tpData(ci).rates)
        rate = tpData(ci).rates(ri).rate;
        Ttp = tpData(ci).rates(ri).Ttp;
        % Throughput summary blocks are normalized into columns like RAW_BER/K/AvgIter.
        % Use the same selector as progress plots.
        [x, y] = selectIterXY(Ttp, model);
        if isempty(x)
            continue;
        end
        [x, y] = sortByX(x, y);

        if strcmpi(model, "AWGN")
            semilogx(ax, x, y, "-o", "LineWidth", 1.2, "MarkerSize", 4);
        else
            plot(ax, x, y, "-o", "LineWidth", 1.2, "MarkerSize", 4);
        end
        leg(end+1) = rate; %#ok<AGROW>
    end

    xlabel(ax, axisName, "Interpreter", "none");
    ylabel(ax, "AvgIter (from Throughput\_summary)", "Interpreter", "none");
    title(ax, figPrefix + cfg + " | Throughput\_summary", "Interpreter", "none");
    if strcmpi(model, "ERR_INJ") && strcmpi(opt.KDirection, "decreasing")
        set(ax, "XDir", "reverse");
    end
    if ~isempty(leg)
        legend(ax, leg, "Interpreter", "none", "Location", "best");
    end
end

%% Plot 4: RBER - K - Throughput (3D, user-defined placeholder)
% This plot is mainly meaningful for ERR_INJ where both RAW_BER and K exist.
% You can define throughput in code by editing this section (kept as placeholder).
throughputFn = @(T) NaN; %#ok<NASGU>
for ci = 1:numel(ferData)
    cfg = ferData(ci).config;
    if isempty(ferData(ci).rates)
        continue;
    end
    [model, ~] = inferModelAndAxis(opt.Model, ferData(ci).rates(1).T);
    if ~strcmpi(model, "ERR_INJ")
        continue;
    end

    f = figure("Name", figPrefix + "Throughput3D_" + cfg, "NumberTitle", "off", "Color", "w");
    ax = axes(f);
    hold(ax, "on");
    grid(ax, "on");

    for ri = 1:numel(ferData(ci).rates)
        rate = ferData(ci).rates(ri).rate;
        T = ferData(ci).rates(ri).T;
        if ~all(ismember(["RAW_BER","K"], string(T.Properties.VariableNames)))
            continue;
        end
        thr = arrayfun(@(~) NaN, (1:height(T))');
        % To enable, replace with:
        % thr = arrayfun(@(i) throughputFn(T(i,:)), (1:height(T))');

        mask = isfinite(T.RAW_BER) & isfinite(T.K) & isfinite(thr);
        if ~any(mask)
            continue;
        end
        scatter3(ax, T.RAW_BER(mask), T.K(mask), thr(mask), 18, "filled", "DisplayName", rate);
    end

    set(ax, "XScale", "log");
    xlabel(ax, "RAW\_BER (RBER)", "Interpreter", "none");
    ylabel(ax, "K (Fixed Errors)", "Interpreter", "none");
    zlabel(ax, "Throughput (user-defined)", "Interpreter", "none");
    title(ax, cfg + " | ERR_INJ | Throughput (placeholder)", "Interpreter", "none");
    legend(ax, "Interpreter", "none", "Location", "best");
end

end

%% Helper Functions
function xlsxPath = resolveSummaryXlsx(p, pattern)
p = string(p);
if strlength(p) == 0
    xlsxPath = "";
    return;
end
if isfile(p)
    xlsxPath = p;
    return;
end
if isfolder(p)
    files = dir(fullfile(p, pattern));
    if isempty(files)
        xlsxPath = "";
        return;
    end
    [~, idx] = max([files.datenum]);
    xlsxPath = string(fullfile(files(idx).folder, files(idx).name));
    return;
end
xlsxPath = "";
end

function [rates, blocksByRate] = loadSummaryWorkbookBlocks(xlsxPath)
if strlength(xlsxPath) == 0 || ~isfile(xlsxPath)
    rates = strings(0, 1);
    blocksByRate = struct();
    return;
end
sheets = getSheetNames(xlsxPath);
rates = string(sheets);
rates = rates(endsWith(rates, "B"));
blocksByRate = struct();
for i = 1:numel(rates)
    rate = rates(i);
    blocks = parseRateSheetBlocks(xlsxPath, rate);
    blocksByRate.(matlab.lang.makeValidName(rate)) = blocks;
end
end

function sheets = getSheetNames(xlsxPath)
try
    sheets = sheetnames(xlsxPath);
catch
    [~, sheets] = xlsfinfo(xlsxPath);
end
end

function blocks = parseRateSheetBlocks(xlsxPath, rateSheet)
C = readcell(xlsxPath, "Sheet", rateSheet);
tmpl = newBlockTemplate();
blocks = tmpl([]);
cur = tmpl;
inBlock = false;

for r = 1:size(C, 1)
    v = C{r, 1};
    if isstring(v) || ischar(v)
        s = string(v);
        if ismissing(s)
            s = "";
        end
    else
        s = "";
    end

    if strlength(s) > 0 && contains(s, "Code Rate:", "IgnoreCase", true)
        if inBlock
            blocks(end+1) = finalizeBlock(cur); %#ok<AGROW>
        end
        cur = initBlockFromHeader(s);
        inBlock = true;
        continue;
    end

    if ~inBlock
        continue;
    end

    % Block terminator: empty metric cell
    if isBlankMetricCell(v)
        blocks(end+1) = finalizeBlock(cur); %#ok<AGROW>
        cur = tmpl;
        inBlock = false;
        continue;
    end

    if ~(isstring(v) || ischar(v))
        % Skip unexpected row shapes inside a block.
        continue;
    end

    metricName = string(v);
    if ismissing(metricName) || strlength(strtrim(metricName)) == 0
        blocks(end+1) = finalizeBlock(cur); %#ok<AGROW>
        cur = tmpl;
        inBlock = false;
        continue;
    end
    rowVals = C(r, 2:end);
    cur.metrics.(matlab.lang.makeValidName(metricName)) = cellRowToNumeric(rowVals);
    cur.metricNames(end+1) = metricName; %#ok<AGROW>
end

if inBlock
    blocks(end+1) = finalizeBlock(cur); %#ok<AGROW>
end
end

function blk = initBlockFromHeader(header)
blk = newBlockTemplate();
blk.header = string(header);
blk.header = replace(blk.header, "｜", "|");

% Parse: "Code Rate: 1952B | Decoder: BF | Channel: AWGN"
pat = "Code\\s*Rate:\\s*(?<rate>[^|]+)\\|\\s*Decoder:\\s*(?<decoder>[^|]+)\\|\\s*Channel:\\s*(?<channel>.+)$";
m = regexp(blk.header, pat, "names", "once", "ignorecase");
if isempty(m)
    [blk.rate, blk.decoder, blk.channel] = parseHeaderFallback(blk.header);
else
    blk.rate = strtrim(string(m.rate));
    blk.decoder = strtrim(string(m.decoder));
    blk.channel = strtrim(string(m.channel));
end
blk.config = canonicalConfigId(blk.decoder, blk.channel);
end

function [rate, decoder, channel] = parseHeaderFallback(header)
rate = "";
decoder = "";
channel = "";
parts = split(string(header), "|");
for i = 1:numel(parts)
    part = strtrim(parts(i));
    if startsWith(part, "Code Rate:", "IgnoreCase", true)
        rate = strtrim(extractAfter(part, ":"));
    elseif startsWith(part, "Decoder:", "IgnoreCase", true)
        decoder = strtrim(extractAfter(part, ":"));
    elseif startsWith(part, "Channel:", "IgnoreCase", true)
        channel = strtrim(extractAfter(part, ":"));
    end
end
end

function cfg = canonicalConfigId(decoder, channel)
dec = strtrim(string(decoder));
dec = regexprep(dec, "\\s+", "_");
ch = upper(strtrim(string(channel)));
if contains(ch, "ERR", "IgnoreCase", true)
    ch = "ERRINJ";
else
    ch = "AWGN";
end
cfg = dec + "_" + ch;
end

function blk = finalizeBlock(blk)
% Normalize metrics into a MATLAB table T.
blk.T = blockMetricsToTable(blk.metrics, blk.metricNames);
end

function x = cellRowToNumeric(rowVals)
rowVals = rowVals(:)';
x = nan(1, numel(rowVals));
for i = 1:numel(rowVals)
    v = rowVals{i};
    if isnumeric(v)
        x(i) = double(v);
    elseif islogical(v)
        x(i) = double(v);
    elseif isstring(v) || ischar(v)
        x(i) = str2double(string(v));
    else
        x(i) = nan;
    end
end
% Trim trailing NaNs (Excel often pads columns).
last = find(isfinite(x), 1, "last");
if isempty(last)
    x = [];
else
    x = x(1:last);
end
end

function T = blockMetricsToTable(metrics, metricNames)
% Build a point-wise table from metric-row vectors.
keys = string(fieldnames(metrics));
if isempty(keys)
    T = table();
    return;
end

% Map normalized metric-name -> original fieldname
normMap = containers.Map("KeyType", "char", "ValueType", "char");
for i = 1:numel(metricNames)
    k = metricNames(i);
    field = matlab.lang.makeValidName(k);
    nk = normMetricName(k);
    if ~isKey(normMap, nk)
        normMap(nk) = field;
    end
end

% Helper to fetch by normalized aliases
getv = @(aliases) pickMetric(metrics, normMap, aliases);

snr = getv(["SNR"]);
kfix = getv(["K"]);
raw = getv(["RAW_BER","RBER","RAWBER"]);
fer = getv(["LDPC_FER","FER"]);
fcw = getv(["FAIL_CW","FAILCW"]);
avg = getv(["AVGITER","AVER_ITER","AVERITER","RETRY_DECODER_AVERAGE_ITERATIONS","RETRY_DECODER_AVERAGE_ITERATION"]);

% Determine length
len = max([numel(snr), numel(kfix), numel(raw), numel(fer), numel(fcw), numel(avg), 0]);
pad = @(v) padToLen(v, len);
snr = pad(snr); kfix = pad(kfix); raw = pad(raw); fer = pad(fer); fcw = pad(fcw); avg = pad(avg);

T = table();
if any(isfinite(snr)); T.SNR = snr(:); end
if any(isfinite(kfix)); T.K = kfix(:); end
if any(isfinite(raw)); T.RAW_BER = raw(:); end
if any(isfinite(fer)); T.LDPC_FER = fer(:); end
if any(isfinite(fcw)); T.FAIL_CW = fcw(:); end
if any(isfinite(avg)); T.AvgIter = avg(:); end
end

function nk = normMetricName(name)
s = upper(string(name));
s = regexprep(s, "[^A-Z0-9]+", "_");
s = regexprep(s, "^_+|_+$", "");
nk = char(s);
end

function v = pickMetric(metrics, normMap, aliases)
v = [];
for a = aliases
    k = normMetricName(a);
    if isKey(normMap, k)
        field = normMap(k);
        v = metrics.(field);
        return;
    end
end
end

function v = padToLen(v, len)
if isempty(v)
    v = nan(len, 1);
    return;
end
v = v(:);
if numel(v) < len
    v(end+1:len, 1) = nan;
elseif numel(v) > len
    v = v(1:len, 1);
end
end

function [rates, cfgs, ferData] = buildFerDataFromBlocks(ratesAll, blocksByRate, ratesIn, cfgsIn, minFailCw, modelOpt)
rates = string(ratesIn);
if isempty(rates)
    rates = ratesAll;
else
    rates = intersect(ratesAll, rates, "stable");
end

% Collect all configs
cfgSet = strings(0, 1);
for i = 1:numel(ratesAll)
    rate = ratesAll(i);
    blocks = blocksByRate.(matlab.lang.makeValidName(rate));
    for b = 1:numel(blocks)
        cfgSet(end+1) = blocks(b).config; %#ok<AGROW>
    end
end
cfgSet = unique(cfgSet, "stable");

cfgs = string(cfgsIn);
if isempty(cfgs)
    cfgs = cfgSet;
else
    % Allow substring match
    keep = false(size(cfgSet));
    for i = 1:numel(cfgSet)
        for j = 1:numel(cfgs)
            if contains(cfgSet(i), string(cfgs(j)), "IgnoreCase", true)
                keep(i) = true;
                break;
            end
        end
    end
    cfgs = cfgSet(keep);
end

ferData = struct();
for ci = 1:numel(cfgs)
    cfg = cfgs(ci);
    ferData(ci).config = cfg; %#ok<AGROW>
    ferData(ci).rates = struct([]);
    for ri = 1:numel(rates)
        rate = rates(ri);
        blocks = blocksByRate.(matlab.lang.makeValidName(rate));
        blk = findBlock(blocks, cfg);
        if isempty(blk)
            continue;
        end
        T = blk.T;
        T = filterFerTable(T, minFailCw, modelOpt);
        if isempty(T) || height(T) == 0
            continue;
        end
        ferData(ci).rates(end+1).rate = rate; %#ok<AGROW>
        ferData(ci).rates(end).T = T;
    end
end
end

function blk = findBlock(blocks, cfg)
blk = [];
for i = 1:numel(blocks)
    if blocks(i).config == cfg
        blk = blocks(i);
        return;
    end
end
end

function T = filterFerTable(T, minFailCw, modelOpt)
if isempty(T)
    return;
end
vars = string(T.Properties.VariableNames);
need = ["RAW_BER","LDPC_FER"];
for k = 1:numel(need)
    if ~any(vars == need(k))
        T = T([],:);
        return;
    end
end
T = T(isfinite(T.RAW_BER) & isfinite(T.LDPC_FER), :);
if minFailCw > 0 && any(vars == "FAIL_CW")
    T = T(isfinite(T.FAIL_CW) & T.FAIL_CW >= minFailCw, :);
end

% If user forces a model, drop incompatible rows/columns.
modelOpt = upper(string(modelOpt));
if modelOpt == "ERR_INJ"
    if any(vars == "K")
        T = T(isfinite(T.K), :);
    else
        T = T([],:);
    end
elseif modelOpt == "AWGN"
    % ok
end
end

function [rates, cfgs, tpData] = buildTpDataFromBlocks(ratesAll, blocksByRate, ratesIn, cfgsIn)
rates = string(ratesIn);
if isempty(rates)
    rates = ratesAll;
else
    rates = intersect(ratesAll, rates, "stable");
end

cfgSet = strings(0, 1);
for i = 1:numel(ratesAll)
    rate = ratesAll(i);
    blocks = blocksByRate.(matlab.lang.makeValidName(rate));
    for b = 1:numel(blocks)
        cfgSet(end+1) = blocks(b).config; %#ok<AGROW>
    end
end
cfgSet = unique(cfgSet, "stable");

cfgs = string(cfgsIn);
if isempty(cfgs)
    cfgs = cfgSet;
else
    keep = false(size(cfgSet));
    for i = 1:numel(cfgSet)
        for j = 1:numel(cfgs)
            if contains(cfgSet(i), string(cfgs(j)), "IgnoreCase", true)
                keep(i) = true;
                break;
            end
        end
    end
    cfgs = cfgSet(keep);
end

tpData = struct();
for ci = 1:numel(cfgs)
    cfg = cfgs(ci);
    tpData(ci).config = cfg; %#ok<AGROW>
    tpData(ci).rates = struct([]);
    for ri = 1:numel(rates)
        rate = rates(ri);
        blocks = blocksByRate.(matlab.lang.makeValidName(rate));
        blk = findBlock(blocks, cfg);
        if isempty(blk)
            continue;
        end
        Ttp = blk.T;
        if isempty(Ttp) || height(Ttp) == 0
            continue;
        end
        tpData(ci).rates(end+1).rate = rate; %#ok<AGROW>
        tpData(ci).rates(end).Ttp = Ttp;
    end
end
end

function [model, axisName] = inferModelAndAxis(modelOpt, T)
modelOpt = upper(string(modelOpt));
if modelOpt ~= "AUTO"
    model = modelOpt;
else
    vars = string(T.Properties.VariableNames);
    if any(vars == "K")
        model = "ERR_INJ";
    else
        model = "AWGN";
    end
end
axisName = "RAW_BER (RBER)";
if strcmpi(model, "ERR_INJ")
    axisName = "K (Fixed Errors)";
end
end

function [model, axisName] = inferModelFromTp(Ttp)
vars = string(Ttp.Properties.VariableNames);
if any(strcmpi(vars, "k"))
    model = "ERR_INJ";
    axisName = "K (Fixed Errors)";
elseif any(strcmpi(vars, "RBER"))
    model = "AWGN";
    axisName = "RBER";
else
    % Fallback: assume first column is axis
    first = vars(1);
    if contains(lower(first), "k")
        model = "ERR_INJ";
        axisName = "K (Fixed Errors)";
    else
        model = "AWGN";
        axisName = "RBER";
    end
end
end

function [x, y] = selectFerXY(T, model)
if strcmpi(model, "ERR_INJ")
    if ismember("K", string(T.Properties.VariableNames))
        x = T.K;
    else
        x = T{:,1};
    end
else
    x = T.RAW_BER;
end
y = T.LDPC_FER;
mask = isfinite(x) & isfinite(y) & (y > 0);
x = x(mask);
y = y(mask);
end

function [x, y] = selectIterXY(T, model)
if ~ismember("AvgIter", string(T.Properties.VariableNames))
    x = [];
    y = [];
    return;
end
if strcmpi(model, "ERR_INJ")
    if ismember("K", string(T.Properties.VariableNames))
        x = T.K;
    else
        x = T{:,1};
    end
else
    x = T.RAW_BER;
end
y = T.AvgIter;
mask = isfinite(x) & isfinite(y);
x = x(mask);
y = y(mask);
end

function [x, y] = selectTpIterXY(Ttp)
vars = string(Ttp.Properties.VariableNames);
if isempty(vars)
    x = [];
    y = [];
    return;
end
% Prefer RBER over SNR when available (matches "RBER/Fix_Err-aver_iter" requirement)
axisVar = "";
if any(strcmpi(vars, "RBER"))
    axisVar = "RBER";
elseif any(strcmpi(vars, "RAW_BER"))
    axisVar = "RAW_BER";
elseif any(strcmpi(vars, "k"))
    axisVar = vars(find(strcmpi(vars, "k"), 1, "first"));
else
    axisVar = vars(1);
end

iterVar = "";
for c = vars
    if lower(c) == "aver_iter" || lower(c) == "avgiter" || lower(c) == "avg_iter"
        iterVar = c;
        break;
    end
end
if strlength(iterVar) == 0
    iterVar = vars(end);
end

x = Ttp.(axisVar);
y = Ttp.(iterVar);
mask = isfinite(x) & isfinite(y);
x = x(mask);
y = y(mask);
end

function [x, y] = sortByX(x, y)
[x, idx] = sort(x(:));
y = y(idx);
end

function tf = isBlankMetricCell(v)
% readcell() may return empty, NaN, or <missing> for blank spreadsheet cells.
tf = false;
if isempty(v)
    tf = true;
    return;
end
try
    if ismissing(v)
        tf = true;
        return;
    end
catch
    % ismissing may not support some types; ignore
end
if isnumeric(v) && isscalar(v) && isnan(v)
    tf = true;
    return;
end
if isstring(v)
    if ismissing(v) || strlength(strtrim(v)) == 0
        tf = true;
        return;
    end
end
if ischar(v)
    if strlength(strtrim(string(v))) == 0
        tf = true;
        return;
    end
end
end

function blk = newBlockTemplate()
% Fixed field order template to avoid "dissimilar structures" errors.
blk = struct();
blk.header = "";
blk.metrics = struct();
blk.metricNames = strings(0, 1);
blk.rate = "";
blk.decoder = "";
blk.channel = "";
blk.config = "";
blk.T = table();
end
