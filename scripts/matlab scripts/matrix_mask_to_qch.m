function qcH = matrix_mask_to_qch(base_matrix, mask_matrix, varargin)
%MATRIX_MASK_TO_QCH Expand base matrix with mask rules to QC parity-check matrix.
%   qcH = MATRIX_MASK_TO_QCH(base_matrix, mask_matrix, 'TmSize', st, ...)
%   Optional name-value pairs:
%     'TmSize'     : number of top rows belonging to submatrix T (required).
%     'CirSize'    : circulant size; defaults to max(shift)+1.
%     'DropLen'    : drop length (pad_bit); defaults to mode of shifts with mask==1.
%     'OutputFile' : path to write the expanded QC-H (0/1) matrix as text.

parser = inputParser;
parser.addRequired('base_matrix', @(x) isnumeric(x));
parser.addRequired('mask_matrix', @(x) isnumeric(x) && all(size(x) == size(base_matrix)));
parser.addParameter('TmSize', [], @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('CirSize', [], @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('DropLen', [], @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('OutputFile', '', @(s) ischar(s) || isstring(s));
parser.parse(base_matrix, mask_matrix, varargin{:});

tm_sz = parser.Results.TmSize;
cir_sz = parser.Results.CirSize;
drop_len = parser.Results.DropLen;
output_file = char(parser.Results.OutputFile);

if isempty(tm_sz)
    error('matrix_mask_to_qch:TmSizeMissing', '''TmSize'' 参数必须提供。');
end

nz_entries = base_matrix >= 0;
if isempty(cir_sz)
    shifts = base_matrix(nz_entries);
    if isempty(shifts)
        error('matrix_mask_to_qch:NoShifts', '无法从 base_matrix 推断 CirSize。');
    end
    cir_sz = max(shifts) + 1;
end

if isempty(drop_len)
    mask_one_shifts = base_matrix(mask_matrix == 1 & nz_entries);
    if isempty(mask_one_shifts)
        drop_len = cir_sz;
    else
        drop_len = mod(mode(mask_one_shifts), cir_sz);
        if drop_len == 0
            drop_len = cir_sz;
        end
    end
end

mask_len = cir_sz - drop_len;
if mask_len < 0
    error('matrix_mask_to_qch:InvalidLengths', '需满足 drop_len <= cir_sz。');
end

[bm_m, bm_n] = size(base_matrix);
hm_m = bm_m * cir_sz - mask_len;
hm_n = bm_n * cir_sz - mask_len;

row_idx = zeros(nnz(nz_entries) * cir_sz, 1);
col_idx = zeros(nnz(nz_entries) * cir_sz, 1);
ptr = 1;

for i = 1:bm_m
    for j = 1:bm_n
        shift = base_matrix(i, j);
        if shift < 0
            continue;
        end
        mask_flag = mask_matrix(i, j);
        if mask_flag == 1
            span = mask_len;
        elseif mask_flag == 2
            span = drop_len;
        else
            span = cir_sz;
        end
        if span <= 0
            continue;
        end

        if i <= tm_sz
            row_base = (i - 1) * cir_sz;
        else
            row_base = (i - 1) * cir_sz - mask_len;
        end
        col_base = (j - 1) * cir_sz;

        k = (0:span - 1).';
        rows = row_base + k + 1;
        cols = col_base + mod(k + shift, cir_sz) + 1;

        idx_end = ptr + span - 1;
        row_idx(ptr:idx_end) = rows;
        col_idx(ptr:idx_end) = cols;
        ptr = idx_end + 1;
    end
end

row_idx = row_idx(1:ptr - 1);
col_idx = col_idx(1:ptr - 1);
val = ones(size(row_idx));

qcH = sparse(row_idx, col_idx, val, hm_m, hm_n);

if ~isempty(strtrim(output_file))
    fid = fopen(output_file, 'w');
    if fid < 0
        error('matrix_mask_to_qch:WriteFailed', '无法创建输出文件 %s。', output_file);
    end
    dense = full(qcH ~= 0);
    for r = 1:size(dense, 1)
        fprintf(fid, '%d ', dense(r, 1:end-1));
        fprintf(fid, '%d\n', dense(r, end));
    end
    fclose(fid);
end
end
