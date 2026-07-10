function plot_autofer_rber_fer(varargin)
%PLOT_AUTOFER_RBER_FER Plot RBER-FER comparison curves from AutoFER xlsx outputs.
%
% Command-style usage:
%   plot_autofer_rber_fer('case_a', 'case_b', 'case_c', ...
%                         '--out-dir', 'plots', ...
%                         '--min-fail-cw', '10');
%
% Name-value usage:
%   plot_autofer_rber_fer({'case_a', 'case_b'}, ...
%                         'out_dir', 'plots', ...
%                         'min_fail_cw', 10);

opts = parse_args(varargin{:});

if isempty(opts.curve_dirs)
    error('plot_autofer_rber_fer:NoCurveDirs', 'At least one curve directory is required.');
end

ensure_dir(opts.out_dir);

curves = struct('dir', {}, 'xlsx', {}, 'label', {}, 'x', {}, 'y', {});
for i = 1:numel(opts.curve_dirs)
    curve_dir = char(opts.curve_dirs(i));
    if ~isfolder(curve_dir)
        error('plot_autofer_rber_fer:BadCurveDir', 'Curve directory not found: %s', curve_dir);
    end

    xlsx_path = find_single_xlsx(curve_dir);
    [x, y] = read_curve_from_xlsx(xlsx_path, opts.min_fail_cw);
    if isempty(x)
        warning('plot_autofer_rber_fer:NoPoints', ...
            'No plottable points after filtering: %s', curve_dir);
        continue;
    end

    curves(end + 1).dir = curve_dir; %#ok<AGROW>
    curves(end).xlsx = xlsx_path;
    curves(end).label = directory_label(curve_dir);
    curves(end).x = x;
    curves(end).y = y;
end

if isempty(curves)
    error('plot_autofer_rber_fer:NoCurves', ...
        'No curves have points satisfying min_fail_cw=%d.', opts.min_fail_cw);
end

fig = figure('Name', 'autofer_rber_fer_compare', 'NumberTitle', 'off', 'Color', 'w');
hold on;

colors = make_colors(numel(curves));
for i = 1:numel(curves)
    semilogy(curves(i).x, curves(i).y, '-o', ...
        'Color', colors(i, :), ...
        'LineWidth', 1.2, ...
        'MarkerSize', 4, ...
        'DisplayName', curves(i).label);
end

set(gca, 'YScale', 'log');
hold off;
xlabel('RBER');
ylabel('FER');
grid on;
legend('show', 'Interpreter', 'none', 'Location', 'eastoutside');
title('AutoFER RBER-FER Comparison', 'Interpreter', 'none');

out_base = fullfile(opts.out_dir, 'autofer_rber_fer_compare');
save_figure(fig, out_base);

fprintf('plotted %d curve(s) to %s\n', numel(curves), opts.out_dir);
end


function opts = parse_args(varargin)
opts = struct();
opts.curve_dirs = strings(0, 1);
opts.out_dir = "plots_autofer_rber_fer";
opts.min_fail_cw = 0;

args = varargin;
if numel(args) == 1 && iscell(args{1})
    args = args{1};
end

i = 1;
while i <= numel(args)
    current = args{i};
    if iscell(current) || (isstring(current) && numel(current) > 1)
        opts.curve_dirs = append_curve_dirs(opts.curve_dirs, current);
        i = i + 1;
        continue;
    end

    key = string(current);
    key = strip(key);
    if startsWith(key, "--")
        key_norm = lower(strrep(extractAfter(key, 2), "-", "_"));
        switch char(key_norm)
            case 'out_dir'
                [opts.out_dir, i] = read_value(args, i);
            case 'min_fail_cw'
                [v, i] = read_value(args, i);
                opts.min_fail_cw = parse_nonnegative_int(v, 'min_fail_cw');
            case 'curve_dirs'
                [v, i] = read_value(args, i);
                opts.curve_dirs = append_curve_dirs(opts.curve_dirs, v);
            otherwise
                error('plot_autofer_rber_fer:BadArg', 'Unknown argument: %s', args{i});
        end
    elseif is_name_value_key(key)
        key_norm = lower(strrep(key, "-", "_"));
        switch char(key_norm)
            case 'out_dir'
                [opts.out_dir, i] = read_value(args, i);
            case 'min_fail_cw'
                [v, i] = read_value(args, i);
                opts.min_fail_cw = parse_nonnegative_int(v, 'min_fail_cw');
            case 'curve_dirs'
                [v, i] = read_value(args, i);
                opts.curve_dirs = append_curve_dirs(opts.curve_dirs, v);
            otherwise
                error('plot_autofer_rber_fer:BadArg', 'Unknown argument: %s', args{i});
        end
    else
        opts.curve_dirs(end + 1, 1) = key; %#ok<AGROW>
        i = i + 1;
    end
end

opts.out_dir = char(opts.out_dir);
end


function tf = is_name_value_key(key)
tf = any(lower(strrep(key, "-", "_")) == ["out_dir", "min_fail_cw", "curve_dirs"]);
end


function [value, next_i] = read_value(args, i)
if i >= numel(args)
    error('plot_autofer_rber_fer:MissingValue', 'Missing value for argument %s', args{i});
end
value = args{i + 1};
next_i = i + 2;
end


function dirs = append_curve_dirs(dirs, value)
if iscell(value)
    vals = string(value(:));
elseif isstring(value) || ischar(value)
    vals = string(value);
else
    error('plot_autofer_rber_fer:BadCurveDirs', 'curve_dirs must be string, char, or cellstr.');
end
vals = vals(:);
vals = vals(strlength(strip(vals)) > 0);
dirs = [dirs; vals]; %#ok<AGROW>
end


function value = parse_nonnegative_int(raw, name)
value = round(str2double(string(raw)));
if isnan(value) || value < 0
    error('plot_autofer_rber_fer:BadNumericArg', '%s must be a non-negative integer.', name);
end
end


function xlsx_path = find_single_xlsx(curve_dir)
d = dir(fullfile(curve_dir, '*.xlsx'));
d = d(~[d.isdir]);
if isempty(d)
    error('plot_autofer_rber_fer:NoXlsx', 'No .xlsx file found in %s', curve_dir);
end
if numel(d) > 1
    names = string({d.name});
    error('plot_autofer_rber_fer:MultipleXlsx', ...
        'Expected exactly one .xlsx in %s, found: %s', curve_dir, strjoin(names, ', '));
end
xlsx_path = fullfile(d(1).folder, d(1).name);
end


function [x, y] = read_curve_from_xlsx(xlsx_path, min_fail_cw)
T = readtable(xlsx_path, 'FileType', 'spreadsheet', 'VariableNamingRule', 'preserve');
headers = T.Properties.VariableNames;

raw_idx = find_column(headers, ["rawber", "rber"]);
fer_idx = find_column(headers, ["ldpcfer", "fer"]);
fail_idx = find_column(headers, ["failcw", "fail", "failcount"]);

if isempty(raw_idx) || isempty(fer_idx)
    error('plot_autofer_rber_fer:MissingColumns', ...
        'Cannot find RAW_BER/RBER and LDPC_FER/FER columns in %s', xlsx_path);
end
if min_fail_cw > 0 && isempty(fail_idx)
    error('plot_autofer_rber_fer:MissingFailCw', ...
        'min_fail_cw requires a FAIL_CW/FAIL/FAILCOUNT column in %s', xlsx_path);
end

raw = to_double(T.(headers{raw_idx}));
fer = to_double(T.(headers{fer_idx}));
keep = isfinite(raw) & isfinite(fer) & raw > 0 & fer > 0;

if min_fail_cw > 0
    fail_cw = to_double(T.(headers{fail_idx}));
    keep = keep & isfinite(fail_cw) & fail_cw >= min_fail_cw;
end

x = raw(keep);
y = fer(keep);
[x, order] = sort(x);
y = y(order);
end


function idx = find_column(headers, candidates)
idx = [];
for i = 1:numel(headers)
    h = normalize_header(headers{i});
    if any(h == candidates)
        idx = i;
        return;
    end
end
end


function h = normalize_header(x)
h = lower(regexprep(string(x), '[^a-zA-Z0-9]+', ''));
end


function out = to_double(v)
if isnumeric(v) || islogical(v)
    out = double(v);
    return;
end
if iscell(v)
    out = nan(numel(v), 1);
    for i = 1:numel(v)
        out(i) = str2double(string(v{i}));
    end
    return;
end
out = str2double(string(v));
end


function label = directory_label(curve_dir)
[~, name] = fileparts(char(strip(string(curve_dir))));
if isempty(name)
    label = char(strip(string(curve_dir)));
else
    label = name;
end
end


function colors = make_colors(n)
if n <= 0
    colors = zeros(0, 3);
elseif exist('turbo', 'builtin') || exist('turbo', 'file')
    colors = turbo(n);
else
    colors = hsv(n);
end
end


function ensure_dir(path_in)
if ~isfolder(path_in)
    mkdir(path_in);
end
end


function save_figure(fig, base_path)
savefig(fig, [base_path '.fig']);
try
    exportgraphics(fig, [base_path '.png'], 'Resolution', 200);
catch
    print(fig, [base_path '.png'], '-dpng', '-r200');
end
end
