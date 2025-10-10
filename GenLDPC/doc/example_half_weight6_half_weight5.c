/*
 * 示例：生成一半列重为6、一半列重为5的LDPC码
 * 演示行重列重均衡机制
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct WD_pair {
    int weight;  // 权重值
    int num;     // 具有该权重的行/列数量
} WD_pair;

typedef struct WD_vector {
    int n;        // 权重分布的种类数
    WD_pair *wd;  // 权重分布数组
} WD_vector;

void print_distribution(const char* name, WD_vector* dist) {
    printf("\n%s 分布:\n", name);
    printf("----------------------------------------\n");
    int total_ones = 0;
    int total_elements = 0;
    
    for (int i = 0; i < dist->n; i++) {
        printf("  组%d: 权重=%d, 数量=%d\n", 
               i, dist->wd[i].weight, dist->wd[i].num);
        total_ones += dist->wd[i].weight * dist->wd[i].num;
        total_elements += dist->wd[i].num;
    }
    
    printf("----------------------------------------\n");
    printf("  总元素数: %d\n", total_elements);
    printf("  总\"1\"数: %d\n", total_ones);
    printf("  平均权重: %.2f\n", (double)total_ones / total_elements);
}

int main() {
    printf("========================================\n");
    printf("  LDPC码行列重均衡配置示例\n");
    printf("  目标: 一半列重6, 一半列重5\n");
    printf("========================================\n");
    
    // ===== 场景1: 简单示例 (20列, 10行) =====
    printf("\n\n场景1: 简单示例\n");
    printf("================\n");
    
    int base_n = 20;  // 总列数
    int base_m = 10;  // 总行数
    int t_size = 6;   // T矩阵大小（用于编码）
    
    printf("基矩阵大小: %d行 × %d列\n", base_m, base_n);
    printf("T矩阵大小: %d\n", t_size);
    
    // 配置列权重分布
    WD_vector* col_dt = malloc(sizeof(WD_vector));
    col_dt->n = 2;  // 2种权重
    col_dt->wd = malloc(2 * sizeof(WD_pair));
    
    // 一半列重6
    col_dt->wd[0].weight = 6;
    col_dt->wd[0].num = base_n / 2;  // 10列
    
    // 一半列重5
    col_dt->wd[1].weight = 5;
    col_dt->wd[1].num = base_n / 2;  // 10列
    
    print_distribution("列权重", col_dt);
    
    // 计算总"1"数（需要考虑T矩阵部分）
    int total_ones_col = 0;
    for (int i = 0; i < col_dt->n; i++) {
        total_ones_col += col_dt->wd[i].weight * col_dt->wd[i].num;
    }
    
    // T矩阵贡献: t_size行，每行在最后t_size列各有1个"1"
    // 最后t_size列的权重已经在col_dt中计入，需要确保行也匹配
    int t_matrix_ones = t_size;
    
    // 注意：最后t_size列的权重需要被分解为：
    // - T矩阵贡献的部分（每列1个"1"来自T矩阵）
    // - 其他行贡献的部分（每列还有(列重-1)个"1"）
    // 所以总"1"数保持不变，只是分配方式不同
    int ones_to_fill = total_ones_col;
    
    printf("\nT矩阵贡献的\"1\"数: %d\n", t_matrix_ones);
    printf("需要额外填充的\"1\"数: %d\n", ones_to_fill);
    
    // 配置行权重分布（自动均衡）
    WD_vector* row_dt = malloc(sizeof(WD_vector));
    row_dt->n = 3;  // 3种权重
    row_dt->wd = malloc(3 * sizeof(WD_pair));
    
    // 计算行权重
    int effective_rows = base_m - 1;  // 最后一行是T矩阵特殊行
    int avg_row_weight = ones_to_fill / effective_rows;
    int remainder = ones_to_fill % effective_rows;
    
    // 大部分行的权重
    row_dt->wd[0].weight = avg_row_weight;
    row_dt->wd[0].num = effective_rows - remainder;
    
    // 少数行的权重（多1）
    row_dt->wd[1].weight = avg_row_weight + 1;
    row_dt->wd[1].num = remainder;
    
    // T矩阵最后一行（只有1个"1"）
    row_dt->wd[2].weight = 1;
    row_dt->wd[2].num = 1;
    
    print_distribution("行权重", row_dt);
    
    // 验证均衡性
    int total_ones_row = 0;
    for (int i = 0; i < row_dt->n; i++) {
        total_ones_row += row_dt->wd[i].weight * row_dt->wd[i].num;
    }
    
    printf("\n✓ 验证结果:\n");
    printf("  列视角总\"1\"数: %d\n", total_ones_col);
    printf("  行视角总\"1\"数: %d\n", total_ones_row);
    
    if (total_ones_col == total_ones_row) {
        printf("  ✓ 均衡检查通过! 行列权重完美匹配!\n");
    } else {
        printf("  ✗ 均衡检查失败! 差值: %d\n", total_ones_col - total_ones_row);
    }
    
    // ===== 场景2: 5G NR标准尺寸示例 =====
    printf("\n\n========================================\n");
    printf("场景2: 5G NR类似配置\n");
    printf("================\n");
    
    base_n = 68;  // 68列
    base_m = 46;  // 46行
    t_size = 42;  // 42×42的T矩阵
    
    printf("基矩阵大小: %d行 × %d列\n", base_m, base_n);
    printf("T矩阵大小: %d\n", t_size);
    printf("码率: %.3f\n", (double)(base_n - base_m) / base_n);
    
    // 重新配置列权重
    col_dt->wd[0].weight = 6;
    col_dt->wd[0].num = base_n / 2;  // 34列
    
    col_dt->wd[1].weight = 5;
    col_dt->wd[1].num = base_n - base_n / 2;  // 34列（处理奇数）
    
    print_distribution("列权重", col_dt);
    
    // 重新计算总"1"数
    total_ones_col = 0;
    for (int i = 0; i < col_dt->n; i++) {
        total_ones_col += col_dt->wd[i].weight * col_dt->wd[i].num;
    }
    
    t_matrix_ones = t_size;
    ones_to_fill = total_ones_col - t_matrix_ones;
    
    printf("\nT矩阵贡献的\"1\"数: %d\n", t_matrix_ones);
    printf("需要额外填充的\"1\"数: %d\n", ones_to_fill);
    
    // 重新配置行权重
    effective_rows = base_m - 1;
    avg_row_weight = ones_to_fill / effective_rows;
    remainder = ones_to_fill % effective_rows;
    
    row_dt->wd[0].weight = avg_row_weight;
    row_dt->wd[0].num = effective_rows - remainder;
    
    row_dt->wd[1].weight = avg_row_weight + 1;
    row_dt->wd[1].num = remainder;
    
    row_dt->wd[2].weight = 1;
    row_dt->wd[2].num = 1;
    
    print_distribution("行权重", row_dt);
    
    // 验证
    total_ones_row = 0;
    for (int i = 0; i < row_dt->n; i++) {
        total_ones_row += row_dt->wd[i].weight * row_dt->wd[i].num;
    }
    
    printf("\n✓ 验证结果:\n");
    printf("  列视角总\"1\"数: %d\n", total_ones_col);
    printf("  行视角总\"1\"数: %d\n", total_ones_row);
    
    if (total_ones_col == total_ones_row) {
        printf("  ✓ 均衡检查通过! 行列权重完美匹配!\n");
    } else {
        printf("  ✗ 均衡检查失败! 差值: %d\n", total_ones_col - total_ones_row);
    }
    
    // ===== 场景3: 更复杂的分布 =====
    printf("\n\n========================================\n");
    printf("场景3: 复杂权重分布\n");
    printf("================\n");
    printf("30%%列重7 + 40%%列重6 + 30%%列重5\n\n");
    
    base_n = 50;
    base_m = 25;
    t_size = 20;
    
    printf("基矩阵大小: %d行 × %d列\n", base_m, base_n);
    
    // 重新分配内存
    free(col_dt->wd);
    col_dt->n = 3;
    col_dt->wd = malloc(3 * sizeof(WD_pair));
    
    col_dt->wd[0].weight = 7;
    col_dt->wd[0].num = 15;  // 30%
    
    col_dt->wd[1].weight = 6;
    col_dt->wd[1].num = 20;  // 40%
    
    col_dt->wd[2].weight = 5;
    col_dt->wd[2].num = 15;  // 30%
    
    print_distribution("列权重", col_dt);
    
    // 计算行权重
    total_ones_col = 0;
    for (int i = 0; i < col_dt->n; i++) {
        total_ones_col += col_dt->wd[i].weight * col_dt->wd[i].num;
    }
    
    t_matrix_ones = t_size;
    ones_to_fill = total_ones_col - t_matrix_ones;
    
    effective_rows = base_m - 1;
    avg_row_weight = ones_to_fill / effective_rows;
    remainder = ones_to_fill % effective_rows;
    
    row_dt->wd[0].weight = avg_row_weight;
    row_dt->wd[0].num = effective_rows - remainder;
    
    row_dt->wd[1].weight = avg_row_weight + 1;
    row_dt->wd[1].num = remainder;
    
    row_dt->wd[2].weight = 1;
    row_dt->wd[2].num = 1;
    
    print_distribution("行权重", row_dt);
    
    // 最终验证
    total_ones_row = 0;
    for (int i = 0; i < row_dt->n; i++) {
        total_ones_row += row_dt->wd[i].weight * row_dt->wd[i].num;
    }
    
    printf("\n✓ 验证结果:\n");
    printf("  列视角总\"1\"数: %d\n", total_ones_col);
    printf("  行视角总\"1\"数: %d\n", total_ones_row);
    
    if (total_ones_col == total_ones_row) {
        printf("  ✓ 均衡检查通过! 行列权重完美匹配!\n");
    } else {
        printf("  ✗ 均衡检查失败! 差值: %d\n", total_ones_col - total_ones_row);
    }
    
    // 清理
    free(col_dt->wd);
    free(col_dt);
    free(row_dt->wd);
    free(row_dt);
    
    printf("\n========================================\n");
    printf("关键要点:\n");
    printf("1. 列权重可以任意指定分布\n");
    printf("2. 行权重通过除法自动计算均衡\n");
    printf("3. 必须满足: Σ(行重) = Σ(列重)\n");
    printf("4. T矩阵部分需要预留配额\n");
    printf("========================================\n\n");
    
    return 0;
}
