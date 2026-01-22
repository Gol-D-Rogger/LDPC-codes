function G = buildFerGridFbcHre(Tall, varargin)
%buildFerGridFbcHre Build FER table with rows=FBC, cols=HREbit, log-log interpolation.
%
% Goal:
%   Given a table Tall that contains at least:
%     - HREbit (numeric)
%     - LDPC_FER (numeric)
%     - some FBC column (numeric)
%   build a grid table:
%     first column  : FBC (ascending)
%     other columns : HREbit = 0..70 (step=1 by default), each cell is FER.
%
%   Missing points are filled by:
%     1) per-HRE column interpolation in log-log coordinates
%        (log10(FBC), log10(FER));
%     2) optional cross-HRE interpolation at fixed FBC in semi-log coordinates
%        (HREbit, log10(FER)) to fill intermediate HRE columns when only
%        sparse HRE simulations exist (e.g., 0/5/10/...).
%     3) final NaN filling (per-HRE): for each HRE column, remaining missing
%        FBC points are filled by log-log linear *extrapolation* from the
%        existing (FBC,FER) curve; extrapolated FER is capped to 1 for the
%        high-FBC direction (once reaching 1, it stays at 1). There is no
%        lower cap in the low-FBC direction.
%     4) optional 2-D monotonic projection (step-like surface): enforce
%        nondecreasing FER along FBC and along HRE by alternating isotonic
%        regression (PAV) projections in log10(FER) domain.
%
% Usage:
%   G = buildFerGridFbcHre(Tall)
%   G = buildFerGridFbcHre(Tall, 'FBCField', 'FBC_noHRE')
%   G = buildFerGridFbcHre(Tall, 'FBCGrid', 800:1200)
%   G = buildFerGridFbcHre(Tall, 'HREGrid', 0:5:70)
%
% Name-value options:
%   'FBCField'   : which column to use as FBC x-axis (default: auto)
%   'FERField'   : which column to use as FER (default: 'LDPC_FER')
%   'HREField'   : which column to use as HRE-bit id (default: 'HREbit')
%   'FBCGrid'    : explicit FBC grid (numeric vector). If empty, uses
%                  integer grid floor(min)..ceil(max).
%   'HREGrid'    : explicit HRE grid (numeric vector). If empty, uses unique
%                  0:70 by default.
%   'Extrapolate': true/false. If true, allows log-log extrapolation outside
%                  the observed FBC range (default: false).
%   'InterpolateAcrossHRE': true/false. If true, after per-HRE interpolation,
%                  fill remaining missing cells by interpolating along HREbit
%                  at fixed FBC in log10(FER) (default: true).
%   'FillNaNsWithMinMax': true/false. If true, fill remaining NaNs using
%                  per-HRE log-log extrapolation with FER capped to 1
%                  (default: true).
%   'EnforceMonotonic': true/false. If true, enforce a step-like monotone
%                  surface: for each HRE, FER is nondecreasing with FBC; for
%                  each FBC, FER is nondecreasing with HRE (default: true).
%   'MonotonicIters': number of alternating projections (default: 5).
%   'FerFloor': minimum FER used for clamping before log operations
%               (default: 1e-12).
%
% Notes:
%   - log-log interpolation requires FBC>0 and FER>0. Rows with FBC<=0 are kept
%     as-is (no interpolation). FER==0 cannot be used as an anchor in log space.

    p = inputParser;
    p.addParameter('FBCField', '', @(s)ischar(s)||isstring(s));
    p.addParameter('FERField', 'LDPC_FER', @(s)ischar(s)||isstring(s));
    p.addParameter('HREField', 'HREbit', @(s)ischar(s)||isstring(s));
    p.addParameter('FBCGrid', [], @(x)isnumeric(x)&&isvector(x));
    p.addParameter('HREGrid', [], @(x)isnumeric(x)&&isvector(x));
    p.addParameter('Extrapolate', false, @(x)islogical(x)&&isscalar(x));
    p.addParameter('InterpolateAcrossHRE', true, @(x)islogical(x)&&isscalar(x));
    p.addParameter('FillNaNsWithMinMax', true, @(x)islogical(x)&&isscalar(x));
    p.addParameter('EnforceMonotonic', true, @(x)islogical(x)&&isscalar(x));
    p.addParameter('MonotonicIters', 5, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
    p.addParameter('FerFloor', 1e-12, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    p.parse(varargin{:});
    opt = p.Results;

    fbcField = char(opt.FBCField);
    ferField = char(opt.FERField);
    hreField = char(opt.HREField);

    mustHave(Tall, hreField);
    mustHave(Tall, ferField);
    if isempty(fbcField)
        fbcField = pickFbcField(Tall);
    end
    mustHave(Tall, fbcField);

    H = Tall.(hreField);
    X = Tall.(fbcField);
    Y = Tall.(ferField);

    if isempty(opt.HREGrid)
        hreGrid = (0:70);
    else
        hreGrid = unique(opt.HREGrid(:).');
        hreGrid = sort(hreGrid);
    end

    if isempty(opt.FBCGrid)
        xmin = floor(min(X(~isnan(X))));
        xmax = ceil(max(X(~isnan(X))));
        if isempty(xmin) || isempty(xmax) || ~isfinite(xmin) || ~isfinite(xmax)
            error('buildFerGridFbcHre:BadFBC', 'Cannot infer FBC range from %s.', fbcField);
        end
        fbcGrid = xmin:xmax;
    else
        fbcGrid = opt.FBCGrid(:);
        fbcGrid = sort(fbcGrid);
    end

    % Initialize grid with NaNs.
    ferGrid = NaN(numel(fbcGrid), numel(hreGrid));

    % For each HRE column, interpolate from the raw observed (X,Y) pairs directly.
    for hi = 1:numel(hreGrid)
        hk = hreGrid(hi);

        sel = (H == hk) & ~isnan(H) & ~isnan(X) & ~isnan(Y);
        if ~any(sel)
            continue;
        end

        xData = X(sel);
        yData = Y(sel);

        % Only use anchors that are valid in log space.
        ok = (xData > 0) & (yData > 0);
        xData = xData(ok);
        yData = yData(ok);
        if numel(xData) < 2
            continue;
        end

        % Collapse duplicate xData (if multiple logs produce the same FBC) by averaging FER.
        [xU, ~, idx] = unique(xData);
        yU = accumarray(idx, yData, [], @mean);

        % Sort by x for interp1.
        [xU, ord] = sort(xU);
        yU = yU(ord);

        xLog = log10(xU);
        yLog = log10(yU);

        % Query on grid (only for FBC>0).
        qOk = (fbcGrid > 0);
        xq = log10(fbcGrid(qOk));
        if opt.Extrapolate
            yqLog = interp1(xLog, yLog, xq, 'linear', 'extrap');
        else
            yqLog = interp1(xLog, yLog, xq, 'linear', NaN);
        end
        ferGrid(qOk, hi) = 10.^yqLog;

        % If there are exact grid matches in raw data, overwrite them with the observed mean.
        % (This keeps reported points exact even if interpolation has numerical drift.)
        [isExact, loc] = ismember(fbcGrid, xU);
        ferGrid(isExact, hi) = yU(loc(isExact));
    end

    % Optional: fill missing columns across HRE at fixed FBC (semi-log in FER).
    if opt.InterpolateAcrossHRE
        for xi = 1:numel(fbcGrid)
            rowFer = ferGrid(xi, :);

            ok = (~isnan(rowFer)) & (rowFer > 0);
            if nnz(ok) < 2
                continue;
            end

            hData = hreGrid(ok);
            yLog = log10(rowFer(ok));

            miss = isnan(rowFer);
            if ~any(miss)
                continue;
            end

            hq = hreGrid(miss);
            if opt.Extrapolate
                yqLog = interp1(hData, yLog, hq, 'linear', 'extrap');
            else
                yqLog = interp1(hData, yLog, hq, 'linear', NaN);
            end
            rowFer(miss) = 10.^yqLog;
            ferGrid(xi, :) = rowFer;
        end
    end

    % Final NaN handling:
    % - Fill leading/trailing NaNs in each column/row with edge min/max.
    % - If an entire row/col is NaN, map it to global min/max based on position.
    if opt.FillNaNsWithMinMax
        ferGrid = fillNaNsMinMax(ferGrid, fbcGrid, hreGrid, opt.FerFloor);
    end

    % Clamp to valid range to avoid log/plotting issues.
    ferGrid = clampFer(ferGrid, opt.FerFloor);

    % Optional: enforce step-like monotonic surface.
    if opt.EnforceMonotonic
        ferGrid = enforceMonotoneSurface(ferGrid, opt.MonotonicIters, opt.FerFloor);
    end

    % Build output table: first column is FBC.
    G = table(fbcGrid, 'VariableNames', {'FBC'});
    for hi = 1:numel(hreGrid)
        varName = matlab.lang.makeValidName(sprintf('HRE_%dbit', hreGrid(hi)));
        G.(varName) = ferGrid(:, hi);
    end
end

function mustHave(T, field)
    if ~any(strcmp(T.Properties.VariableNames, field))
        error('buildFerGridFbcHre:MissingField', 'Missing required field: %s', field);
    end
end

function fbcField = pickFbcField(T)
% Heuristic priority (project-specific):
%   Prefer FBC_bit whenever present (user-defined FBC axis).
%   Otherwise fall back to other project-specific FBC proxies.
    candidates = { ...
        'FBC_bit', ...
        'FBC_noHRE', ...
        'FBC_TAWGN_x1_x2', ...
        'FBC_total_x0_x2'};
    fbcField = '';
    for i = 1:numel(candidates)
        if any(strcmp(T.Properties.VariableNames, candidates{i}))
            fbcField = candidates{i};
            return;
        end
    end
    error('buildFerGridFbcHre:NoFBCField', ...
        'Cannot auto-pick an FBC field. Please pass ''FBCField'' explicitly.');
end

function A = fillNaNsMinMax(A, fbcGrid, hreGrid, ferFloor)
% Fill remaining NaNs per HRE column by log-log linear extrapolation.
% - High-FBC direction: extrapolate and cap FER to 1 (then stays at 1).
% - Low-FBC direction: extrapolate with no lower cap.
%
% Note: hreGrid is currently unused but kept for interface stability.
    %#ok<NASGU>  % hreGrid intentionally unused
    if nargin < 4 || isempty(ferFloor)
        ferFloor = 1e-12;
    end

    if ~any(isfinite(A(:))) || isempty(fbcGrid)
        return;
    end

    fbcGrid = fbcGrid(:);
    okX = isfinite(fbcGrid) & (fbcGrid > 0);
    if ~any(okX)
        return;
    end

    idxX = find(okX);

    for j = 1:size(A, 2)
        col = A(:, j);
        if numel(col) ~= numel(fbcGrid)
            error('buildFerGridFbcHre:SizeMismatch', ...
                'Size mismatch in fillNaNsMinMax: numel(FBCGrid)=%d but size(FERGrid,1)=%d.', ...
                numel(fbcGrid), numel(col));
        end

        valid = okX & isfinite(col) & (col > 0);
        if ~any(valid)
            continue;
        end

        xLog = log10(fbcGrid(valid));
        yLog = log10(col(valid));
        xLog = xLog(:);
        yLog = yLog(:);

        % Only fill missing/invalid entries; do not overwrite existing valid data.
        need = (~isfinite(col)) | (col <= 0);
        idxFill = find(need & okX);
        if isempty(idxFill)
            continue;
        end

        xqFillLog = log10(fbcGrid(idxFill));

        % Need at least 2 anchors for a slope; with 1 anchor, use constant fill.
        if numel(xLog) == 1
            yqFillLog = repmat(yLog(1), size(xqFillLog));
        else
            [xLog, ord] = sort(xLog);
            yLog = yLog(ord);
            yqFillLog = interp1(xLog, yLog, xqFillLog, 'linear', 'extrap');
        end

        yqFillLog = yqFillLog(:);
        if numel(yqFillLog) ~= numel(idxFill)
            error('buildFerGridFbcHre:InterpSizeMismatch', ...
                'Interpolation size mismatch: numel(idxFill)=%d but numel(yqFillLog)=%d.', ...
                numel(idxFill), numel(yqFillLog));
        end

        fillVals = 10.^yqFillLog;
        fillVals(fillVals > 1) = 1;
        fillVals(fillVals < ferFloor) = ferFloor;
        col(idxFill) = fillVals;
        A(:, j) = col;
    end
end

function A = clampFer(A, ferFloor)
% Clamp FER to [ferFloor, 1] where finite; keep NaNs as-is.
    if nargin < 2 || isempty(ferFloor)
        ferFloor = 1e-12;
    end
    finite = isfinite(A);
    A(finite & A > 1) = 1;
    A(finite & A < ferFloor) = ferFloor;
end

function A = enforceMonotoneSurface(A, iters, ferFloor)
% Enforce nondecreasing FER down rows (FBC) and across cols (HRE).
%
% Uses alternating projections via isotonic regression (PAV) in log10(FER)
% domain, which typically avoids the "single outlier propagates everywhere"
% behavior of cummax.
    if nargin < 2 || isempty(iters)
        iters = 5;
    end
    if nargin < 3 || isempty(ferFloor)
        ferFloor = 1e-12;
    end
    if iters <= 0
        return;
    end

    % Work in log domain for "decade-wise" smoothness.
    A = clampFer(A, ferFloor);
    Z = log10(A);

    % Remaining NaNs (if any) are left untouched.
    nanMask = isnan(Z);

    for k = 1:iters
        % Column-wise isotonic (along FBC)
        for j = 1:size(Z, 2)
            zcol = Z(:, j);
            ok = ~isnan(zcol);
            if nnz(ok) >= 2
                zcol(ok) = isotonicNondecreasing(zcol(ok));
                Z(:, j) = zcol;
            end
        end

        % Row-wise isotonic (along HRE)
        for i = 1:size(Z, 1)
            zrow = Z(i, :).';
            ok = ~isnan(zrow);
            if nnz(ok) >= 2
                zrow(ok) = isotonicNondecreasing(zrow(ok));
                Z(i, :) = zrow.';
            end
        end
    end

    Z(nanMask) = NaN;
    A = 10.^Z;
    A = clampFer(A, ferFloor);
end

function y = isotonicNondecreasing(y)
% Unweighted isotonic regression (PAV) for nondecreasing sequence.
% Minimizes sum (y - yhat)^2 subject to yhat(1)<=...<=yhat(n).
    y = y(:);
    n = numel(y);
    if n <= 1
        return;
    end

    % PAV blocks
    v = y;                 % block values
    w = ones(n, 1);        % block weights
    start = (1:n).';       % block start index
    stop  = (1:n).';       % block stop index
    m = n;                 % number of blocks

    i = 1;
    while i < m
        if v(i) <= v(i+1)
            i = i + 1;
            continue;
        end
        % merge i and i+1
        newW = w(i) + w(i+1);
        newV = (w(i)*v(i) + w(i+1)*v(i+1)) / newW;
        v(i) = newV;
        w(i) = newW;
        stop(i) = stop(i+1);

        % delete block i+1 by shifting left
        if i+1 < m
            v(i+1:m-1) = v(i+2:m);
            w(i+1:m-1) = w(i+2:m);
            start(i+1:m-1) = start(i+2:m);
            stop(i+1:m-1)  = stop(i+2:m);
        end
        m = m - 1;

        % backtrack if needed
        if i > 1
            i = i - 1;
        end
    end

    % expand blocks
    yhat = zeros(n, 1);
    for b = 1:m
        yhat(start(b):stop(b)) = v(b);
    end
    y = yhat;
end
