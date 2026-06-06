function resultTable = plot_latency_from_throughput_logs(logPath, xAxisType, strobes, userdataB, parityB)
%PLOT_LATENCY_FROM_THROUGHPUT_LOGS Plot latency curves from throughput logs.
%   resultTable = plot_latency_from_throughput_logs(logPath, xAxisType, ...
%       strobes, userdataB, parityB)
%
%   logPath   : directory containing *.log files, or one .log file
%   xAxisType : "RBER" or "RBE"
%   strobes   : strobe count used in legend
%   userdataB : userdata bytes
%   parityB   : parity bytes

% User-editable display label.
iterLabel = '1024,512';

% User-editable latency formula. itr_num is the parsed iteration number.
latencyNs = @(itr_num) (331+2*12)*1000/525*itr_num*2;

xAxisType = upper(string(xAxisType));
if xAxisType ~= "RBER" && xAxisType ~= "RBE"
    error('plot_latency_from_throughput_logs:BadXAxis', 'xAxisType must be RBER or RBE.');
end

logFiles = collect_log_files(logPath);
if isempty(logFiles)
    error('plot_latency_from_throughput_logs:NoLogs', 'No .log files found under %s.', string(logPath));
end

rberVals = [];
rbeVals = [];
itr99Vals = [];
itr9999Vals = [];
itr999999Vals = [];
skipped = strings(0, 1);

for i = 1:numel(logFiles)
    logFile = logFiles(i);
    parsed = parse_one_log(logFile);
    if ~parsed.valid
        skipped(end + 1, 1) = string(logFile); %#ok<AGROW>
        continue;
    end

    rberVals(end + 1, 1) = parsed.rber; %#ok<AGROW>
    rbeVals(end + 1, 1) = round(parsed.rber * (userdataB + parityB) * 8); %#ok<AGROW>
    itr99Vals(end + 1, 1) = parsed.itr99; %#ok<AGROW>
    itr9999Vals(end + 1, 1) = parsed.itr9999; %#ok<AGROW>
    itr999999Vals(end + 1, 1) = parsed.itr999999; %#ok<AGROW>
end

if isempty(rberVals)
    error('plot_latency_from_throughput_logs:NoValidLogs', 'No valid logs found. All logs were skipped.');
end

if xAxisType == "RBER"
    xVals = rberVals;
    xLabelText = 'RBER';
else
    xVals = rbeVals;
    xLabelText = 'RBE';
end

lat99 = latencyNs(itr99Vals);
lat9999 = latencyNs(itr9999Vals);
lat999999 = latencyNs(itr999999Vals);

[xVals, order] = sort(xVals);
rberVals = rberVals(order);
rbeVals = rbeVals(order);
itr99Vals = itr99Vals(order);
itr9999Vals = itr9999Vals(order);
itr999999Vals = itr999999Vals(order);
lat99 = lat99(order);
lat9999 = lat9999(order);
lat999999 = lat999999(order);

figure;
plot(xVals, lat99, '-o', 'LineWidth', 1.8, 'MarkerSize', 6);
hold on;
plot(xVals, lat9999, '-s', 'LineWidth', 1.8, 'MarkerSize', 6);
plot(xVals, lat999999, '-^', 'LineWidth', 1.8, 'MarkerSize', 6);
hold off;

grid on;
xlabel(xLabelText);
ylabel('latency(ns)');
title(sprintf('Latency vs %s', xLabelText));
legend( ...
    sprintf('%d strobes, Iter: %s, Userdata %dB, Parity %dB, 99%%', strobes, iterLabel, userdataB, parityB), ...
    sprintf('%d strobes, Iter: %s, Userdata %dB, Parity %dB, 99.99%%', strobes, iterLabel, userdataB, parityB), ...
    sprintf('%d strobes, Iter: %s, Userdata %dB, Parity %dB, 99.9999%%', strobes, iterLabel, userdataB, parityB), ...
    'Location', 'best');

resultTable = table( ...
    rberVals, rbeVals, itr99Vals, itr9999Vals, itr999999Vals, lat99, lat9999, lat999999, ...
    'VariableNames', {'RBER', 'RBE', 'Iter99', 'Iter9999', 'Iter999999', 'Latency99_ns', 'Latency9999_ns', 'Latency999999_ns'});

outDir = output_dir_for(logPath);
outPng = fullfile(outDir, sprintf('latency_vs_%s.png', lower(char(xAxisType))));
outCsv = fullfile(outDir, sprintf('latency_vs_%s.csv', lower(char(xAxisType))));
saveas(gcf, outPng);
writetable(resultTable, outCsv);

fprintf('Parsed %d valid log(s). Skipped %d incomplete/invalid log(s).\n', height(resultTable), numel(skipped));
fprintf('Figure saved: %s\n', outPng);
fprintf('Table saved : %s\n', outCsv);
if ~isempty(skipped)
    fprintf('Skipped logs:\n');
    for i = 1:numel(skipped)
        fprintf('  %s\n', skipped(i));
    end
end
end


function logFiles = collect_log_files(logPath)
p = string(logPath);
if isfile(p)
    logFiles = string(p);
    return;
end
if ~isfolder(p)
    error('plot_latency_from_throughput_logs:BadPath', 'Path does not exist: %s', p);
end

d = dir(fullfile(p, '*.log'));
logFiles = strings(numel(d), 1);
for i = 1:numel(d)
    logFiles(i) = string(fullfile(d(i).folder, d(i).name));
end
end


function outDir = output_dir_for(logPath)
p = string(logPath);
if isfile(p)
    outDir = char(fileparts(p));
else
    outDir = char(p);
end
end


function parsed = parse_one_log(logFile)
parsed = struct('valid', false, 'rber', NaN, 'itr99', NaN, 'itr9999', NaN, 'itr999999', NaN);
try
    txt = fileread(logFile);
catch
    return;
end

rber = last_number(txt, '\[STATISTICS\].*RAW\s+BER\s*:\s*([0-9.eE+-]+)');
itr99 = last_number(txt, '\[STATISTICS\].*iteration above\s+99%\s*:?\s*([0-9.eE+-]+)');
itr9999 = last_number(txt, '\[STATISTICS\].*iteration above\s+99\.99%\s*:?\s*([0-9.eE+-]+)');
itr999999 = last_number(txt, '\[STATISTICS\].*iteration above\s+99\.9999%\s*:?\s*([0-9.eE+-]+)');

if any(isnan([rber, itr99, itr9999, itr999999]))
    return;
end

parsed.valid = true;
parsed.rber = rber;
parsed.itr99 = itr99;
parsed.itr9999 = itr9999;
parsed.itr999999 = itr999999;
end


function value = last_number(txt, pattern)
tokens = regexp(txt, pattern, 'tokens');
if isempty(tokens)
    value = NaN;
    return;
end
value = str2double(tokens{end}{1});
end
