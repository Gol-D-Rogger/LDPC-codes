function qcH = matrix_mask_to_qch(base_matrix, mask_matrix, output_file)
%MATRIX_MASK_TO_QCH 将 mask 为 1 的位置在 base_matrix 中置为 -1。
%   qcH = MATRIX_MASK_TO_QCH(base_matrix, mask_matrix, output_file)
%   会生成更新后的矩阵并写入 output_file，输出文本中每个元素占 4 位。

qcH = base_matrix;
qcH(mask_matrix == 1) = -1;

fid = fopen(output_file, 'w');
if fid < 0
    error('matrix_mask_to_qch:WriteFailed', '无法创建输出文件 %s。', output_file);
end

[m, n] = size(qcH);
for r = 1:m
    fprintf(fid, '%4d', qcH(r, :));
    fprintf(fid, '\n');
end

fclose(fid);
end
