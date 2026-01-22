function [Tall, fig] = plotRberFerAllHre(rootDir, prefix, mode, varargin)
%plotRberFerAllHre Plot RBER-FER (LDPC_FER) curves for all HREbit folders.
%
% This script reads logs under subfolders like:
%   rootDir/0bit, rootDir/1bit, ..., rootDir/70bit
% using extractStatsFromHreFolders(), and plots one curve per HREbit.
%
% Usage:
%   plotRberFerAllHre(rootDir, prefix, mode)
%   plotRberFerAllHre(rootDir, prefix, mode, 'RBERField','TheoRBER')
%   plotRberFerAllHre(rootDir, prefix, mode, 'SavePath','fig_rber_fer.png')
%
% Inputs:
%   rootDir : parent directory containing "<number>bit" subfolders
%   prefix  : log filename prefix (same meaning as extractStatsFromLogs)
%   mode    : log parsing mode (1=extended, 2=minimal)
%
% Name-value options:
%   'RBERField' : 'TheoRBER' (default) or 'RAW_BER'
%   'SavePath' : if non-empty, saves the figure (format inferred from extension)
%
% Outputs:
%   Tall : concatenated table from extractStatsFromHreFolders()
%   fig  : handle to the created figure

    if nargin < 1 || isempty(rootDir), rootDir = '.'; end
    if nargin < 2, prefix = ''; end
    if nargin < 3 || isempty(mode), mode = 1; end

    p = inputParser;
    p.addParameter('RBERField', 'TheoRBER', @(s)ischar(s)||isstring(s));
    p.addParameter('SavePath', '', @(s)ischar(s)||isstring(s));
    p.parse(varargin{:});
    opt = p.Results;

    rberField = char(opt.RBERField);

    Tall = extractStatsFromHreFolders(rootDir, prefix, mode);
    if ~any(strcmp(Tall.Properties.VariableNames, rberField))
        error('plotRberFerAllHre:BadRBERField', 'Tall does not contain RBERField=%s.', rberField);
    end
    if ~any(strcmp(Tall.Properties.VariableNames, 'LDPC_FER'))
        error('plotRberFerAllHre:MissingFER', 'Tall does not contain LDPC_FER.');
    end

    % Keep only valid points for plotting in log-log.
    ok = isfinite(Tall.(rberField)) & isfinite(Tall.LDPC_FER) & ...
         (Tall.(rberField) > 0) & (Tall.LDPC_FER > 0);
    Tall = Tall(ok, :);

    hreList = unique(Tall.HREbit);
    hreList = sort(hreList(:));
    if isempty(hreList)
        error('plotRberFerAllHre:NoData', 'No valid (RBER>0, FER>0) points found to plot.');
    end

    fig = figure('Color','w');
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on');

    cmap = localColormap(max(numel(hreList), 2));

    for i = 1:numel(hreList)
        hre = hreList(i);
        sel = Tall.HREbit == hre;
        x = Tall.(rberField)(sel);
        y = Tall.LDPC_FER(sel);

        % Aggregate duplicate x by averaging y (robust to repeated runs).
        [xU, ~, idx] = unique(x);
        yU = accumarray(idx, y, [], @mean);
        [xU, ord] = sort(xU);
        yU = yU(ord);

        c = cmap(min(i, size(cmap,1)), :);
        loglog(ax, xU, yU, '-', 'Color', c, 'LineWidth', 0.9);
    end

    % Optional baseline curve overlay (user-pasted arrays).
    [baseEnable, baseX, baseY, baseLabel, baseStyle] = userBaseline();
    if baseEnable
        baseX = baseX(:);
        baseY = baseY(:);
        okb = isfinite(baseX) & isfinite(baseY) & (baseX > 0) & (baseY > 0);
        baseX = baseX(okb);
        baseY = baseY(okb);
        if numel(baseX) < 2
            error('plotRberFerAllHre:BadBaseline', ...
                'Baseline enabled but has <2 valid points. Please edit userBaseline() in this file.');
        end
        [xU, ~, idx] = unique(baseX);
        yU = accumarray(idx, baseY, [], @mean);
        [xU, ord] = sort(xU);
        yU = yU(ord);
        loglog(ax, xU, yU, baseStyle, 'LineWidth', 2.0, 'DisplayName', baseLabel);
        legend(ax, 'show', 'Location', 'southwest');
    end

    grid(ax, 'on');
    ax.XScale = 'log';
    ax.YScale = 'log';
    xlabel(ax, sprintf('%s (bit error rate)', rberField), 'Interpreter','none');
    ylabel(ax, 'LDPC\_FER', 'Interpreter','none');
    title(ax, 'RBER-FER Curves Across HRE Bits', 'Interpreter','none');

    % Use colorbar as HREbit index proxy (legend would be too large).
    colormap(ax, cmap);
    cb = colorbar(ax);
    cb.Label.String = 'HREbit';
    cb.Ticks = linspace(0, 1, min(8, numel(hreList)));
    cb.TickLabels = arrayfun(@(v) sprintf('%d', v), ...
        round(linspace(hreList(1), hreList(end), numel(cb.Ticks))), ...
        'UniformOutput', false);

    hold(ax, 'off');

    savePath = char(opt.SavePath);
    if ~isempty(savePath)
        exportgraphics(fig, savePath);
    end
end

function cmap = localColormap(n)
% Choose a visually distinct colormap, with fallback for older MATLAB versions.
    if exist('turbo', 'file') == 2
        cmap = turbo(n);
    else
        cmap = parula(n);
    end
end

function [enable, x, y, label, style] = userBaseline()
% USER BASELINE MODULE
%
% If you want to overlay a baseline curve, edit this function:
%   - set enable=true
%   - paste baseline RBER into x
%   - paste baseline FER into y
%
% IMPORTANT:
%   - x must use the same meaning as the plot x-axis (RBERField). This script
%     does not convert SNR->RBER; it only plots (x,y) directly.
%   - x and y must be positive for log-log plotting.

    enable = false;
    label = "Baseline";
    style = 'k--';

    % Example (replace with your data):
    % x = [1e-2 5e-3 2e-3 1e-3];
    % y = [1e-1 1e-2 1e-4 1e-6];
    x = [];
    y = [];
end
