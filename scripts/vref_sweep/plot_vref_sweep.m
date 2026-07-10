function plot_vref_sweep(varargin)
%PLOT_VREF_SWEEP Plot VREF sweep RBER-FER curves from AutoFER progress xlsx.
%
% Command-style usage:
%   plot_vref_sweep('--input-root', 'output/vref_sweep', ...
%                   '--out-dir', 'output/vref_sweep/plots_matlab', ...
%                   '--min-fail-cw', 10)
%
% Name-value usage:
%   plot_vref_sweep('input_root', 'output/vref_sweep', ...
%                   'min_fail_cw', 10)

opts = parse_args(varargin{:});

input_root = char(opts.input_root);
if ~isfolder(input_root)
    error('plot_vref_sweep:InputRootNotFound', 'input root not found: %s', input_root);
end

if strlength(opts.metadata) > 0
    metadata_path = char(opts.metadata);
else
    metadata_path = fullfile(input_root, 'cases.csv');
end

if strlength(opts.out_dir) > 0
    out_dir = char(opts.out_dir);
else
    out_dir = fullfile(input_root, 'plots_matlab');
end

metadata = load_metadata(metadata_path);
xlsx_files = discover_xlsx(input_root);
if isempty(xlsx_files)
    error('plot_vref_sweep:NoXlsx', 'No *_progress.xlsx files found under %s', input_root);
end

curves = struct('path', {}, 'label', {}, 'x', {}, 'y', {}, 'case_name', {});
for i = 1:numel(xlsx_files)
    xlsx_path = xlsx_files{i};
    [x, y] = read_curve_from_xlsx(xlsx_path, opts.complete_only, opts.min_fail_cw);
    if isempty(x)
        continue;
    end

    cname = case_name_from_xlsx(xlsx_path);
    lbl = label_for(xlsx_path, cname, metadata);
    curves(end + 1).path = xlsx_path; %#ok<AGROW>
    curves(end).label = lbl;
    curves(end).x = x;
    curves(end).y = y;
    curves(end).case_name = cname;

end

if isempty(curves)
    error('plot_vref_sweep:NoPoints', 'No plottable RBER/FER points found');
end

ensure_dir(out_dir);
fig = figure('Name', 'all_vref_overlay', 'NumberTitle', 'off', 'Color', 'w');
hold on;

colors = make_colors(numel(curves));
for i = 1:numel(curves)
    if is_baseline_curve(curves(i))
        line_style = '--o';
    else
        line_style = '-o';
    end
    semilogy(curves(i).x, curves(i).y, line_style, ...
        'Color', colors(i, :), 'LineWidth', 1.1, 'MarkerSize', 4, ...
        'DisplayName', curves(i).label);
end
set(gca, 'YScale', 'log');
hold off;
xlabel('RBER');
ylabel('FER');
grid on;
legend('show', 'Interpreter', 'none', 'Location', 'eastoutside');
title('VREF Sweep RBER-FER', 'Interpreter', 'none');
save_figure(fig, fullfile(out_dir, 'all_vref_overlay'));

fprintf('plotted %d curve(s) to %s\n', numel(curves), out_dir);
end


function tf = is_baseline_curve(curve)
cname = lower(string(curve.case_name));
lbl = lower(string(curve.label));
tf = cname == "baseline" || startsWith(lbl, "baseline");
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


function opts = parse_args(varargin)
opts = struct();
opts.input_root = "output/vref_sweep";
opts.metadata = "";
opts.out_dir = "";
opts.complete_only = false;
opts.min_fail_cw = 0;

args = varargin;
if numel(args) == 1 && iscell(args{1})
    args = args{1};
end

i = 1;
while i <= numel(args)
    key = string(args{i});
    key = strip(key);
    if startsWith(key, "--")
        key = extractAfter(key, 2);
    end
    key = strrep(key, "-", "_");

    switch lower(char(key))
        case 'input_root'
            [opts.input_root, i] = read_value(args, i);
        case 'metadata'
            [opts.metadata, i] = read_value(args, i);
        case 'out_dir'
            [opts.out_dir, i] = read_value(args, i);
        case 'complete_only'
            if i == numel(args) || startsWith(string(args{i + 1}), "--")
                opts.complete_only = true;
                i = i + 1;
            else
                [v, i] = read_value(args, i);
                opts.complete_only = parse_bool(v);
            end
        case 'min_fail_cw'
            [v, i] = read_value(args, i);
            opts.min_fail_cw = max(0, round(str2double(string(v))));
            if isnan(opts.min_fail_cw)
                error('plot_vref_sweep:BadMinFailCw', 'min_fail_cw must be numeric');
            end
        otherwise
            error('plot_vref_sweep:BadArg', 'Unknown argument: %s', args{i});
    end
end
end


function [value, next_i] = read_value(args, i)
if i >= numel(args)
    error('plot_vref_sweep:MissingValue', 'Missing value for argument %s', args{i});
end
value = string(args{i + 1});
next_i = i + 2;
end


function v = parse_bool(x)
s = lower(strip(string(x)));
v = any(s == ["true", "1", "yes", "y", "on"]);
end


function files = discover_xlsx(input_root)
d = dir(fullfile(input_root, '**', '*_progress.xlsx'));
files = cell(numel(d), 1);
for i = 1:numel(d)
    files{i} = fullfile(d(i).folder, d(i).name);
end
files = sort(files);
end


function metadata = load_metadata(metadata_path)
metadata = struct('case_name', {}, 'progress_xlsx', {}, 'vref_label', {}, 'vref_values', {});
if ~isfile(metadata_path)
    return;
end

T = readtable(metadata_path, 'FileType', 'text', 'Delimiter', ',', 'VariableNamingRule', 'preserve');
names = T.Properties.VariableNames;
for i = 1:height(T)
    rec = struct();
    rec.case_name = string(get_table_value(T, names, i, 'case_name'));
    rec.progress_xlsx = string(get_table_value(T, names, i, 'progress_xlsx'));
    rec.vref_label = string(get_table_value(T, names, i, 'vref_label'));
    rec.vref_values = string(get_table_value(T, names, i, 'vref_values'));
    metadata(end + 1) = rec; %#ok<AGROW>
end
end


function value = get_table_value(T, names, row_idx, col_name)
idx = find(strcmp(names, col_name), 1);
if isempty(idx)
    value = "";
    return;
end
col = T.(names{idx});
value = col(row_idx);
if iscell(value)
    value = value{1};
end
if ismissing(value)
    value = "";
end
end


function [x, y] = read_curve_from_xlsx(xlsx_path, complete_only, min_fail_cw)
T = readtable(xlsx_path, 'FileType', 'spreadsheet', 'VariableNamingRule', 'preserve');
headers = T.Properties.VariableNames;

raw_idx = find_column(headers, ["rawber", "rber"]);
fer_idx = find_column(headers, ["ldpcfer", "fer"]);
fail_idx = find_column(headers, ["failcw", "fail", "failcount"]);
complete_idx = find_column(headers, ["iscomplete", "complete"]);

if isempty(raw_idx) || isempty(fer_idx)
    error('plot_vref_sweep:MissingColumns', 'Cannot find RAW_BER/LDPC_FER columns in %s', xlsx_path);
end
if min_fail_cw > 0 && isempty(fail_idx)
    error('plot_vref_sweep:MissingFailCw', '--min-fail-cw requires a FAIL_CW column in %s', xlsx_path);
end

raw = to_double(T.(headers{raw_idx}));
fer = to_double(T.(headers{fer_idx}));
keep = isfinite(raw) & isfinite(fer) & raw > 0 & fer > 0;

if complete_only && ~isempty(complete_idx)
    complete = to_double(T.(headers{complete_idx}));
    keep = keep & complete == 1;
end

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


function cname = case_name_from_xlsx(xlsx_path)
[~, stem, ~] = fileparts(xlsx_path);
suffix = '_progress';
if endsWith(stem, suffix)
    cname = extractBefore(string(stem), strlength(stem) - strlength(suffix) + 1);
else
    cname = string(stem);
end
cname = char(cname);
end


function lbl = label_for(xlsx_path, cname, metadata)
abs_xlsx = string(abs_path(xlsx_path));
for i = 1:numel(metadata)
    if metadata(i).progress_xlsx ~= ""
        if string(abs_path(metadata(i).progress_xlsx)) == abs_xlsx
            lbl = make_label(cname, metadata(i));
            return;
        end
    end
end
for i = 1:numel(metadata)
    if metadata(i).case_name == string(cname)
        lbl = make_label(cname, metadata(i));
        return;
    end
end
lbl = make_vref_label("", cname);
end


function lbl = make_label(cname, rec)
vref = rec.vref_label;
if vref == ""
    vref = rec.vref_values;
end
lbl = make_vref_label(vref, cname);
end


function lbl = make_vref_label(value, fallback)
vref = strip(string(value));
if vref == ""
    vref = string(fallback);
end
if startsWith(lower(vref), "vref=")
    lbl = char(vref);
else
    lbl = char("VREF=" + vref);
end
end


function p = abs_path(path_in)
p = char(path_in);
if exist('javaObject', 'builtin') || exist('javaObject', 'file')
    try
        p = char(java.io.File(p).getCanonicalPath());
        return;
    catch
    end
end
if ~isfolder(p) && ~isfile(p)
    [folder, name, ext] = fileparts(p);
    if isempty(folder)
        folder = pwd;
    end
    p = fullfile(folder, [name ext]);
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
