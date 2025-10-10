/*
 * 行列重均衡机制详细演示
 * 重点：展示如何生成一半列重6、一半列重5的LDPC码
 */

#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int weight;  // 权重值
    int num;     // 该权重的数量
} WD_pair;

typedef struct {
    int n;        // 分布种类数
    WD_pair *wd;  // 权重对数组
} WD_vector;

void print_detail(const char* title, WD_vector* dist, int total_elements) {
    printf("\n%s:\n", title);
    printf("--------------------------------------------------\n");
    int total_ones = 0;
    
    for (int i = 0; i < dist->n; i++) {
        printf("  第%d组: 权重=%2d × 数量=%3d = %4d个\"1\"\n", 
               i+1, dist->wd[i].weight, dist->wd[i].num,
               dist->wd[i].weight * dist->wd[i].num);
        total_ones += dist->wd[i].weight * dist->wd[i].num;
    }
    
    printf("--------------------------------------------------\n");
    printf("  合计: %d个元素, %d个\"1\", 平均权重=%.2f\n", 
           total_elements, total_ones, (double)total_ones/total_elements);
}

int main() {
    printf("\n");
    printf("╔═══════════════════════════════════════════════════╗\n");
    printf("║     LDPC码行列重均衡机制完整演示                 ║\n");
    printf("║     目标: 一半列重6 + 一半列重5                  ║\n");
    printf("╚═══════════════════════════════════════════════════╝\n");
    
    // ============= 配置参数 =============
    int base_n = 40;  // 基矩阵列数
    int base_m = 20;  // 基矩阵行数
    int t_size = 10;  // T矩阵大小（用于系统编码）
    
    printf("\n【第1步】 基本参数设置\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("  • 基矩阵尺寸: %d行 × %d列\n", base_m, base_n);
    printf("  • T矩阵尺寸: %d × %d (用于系统编码的双对角结构)\n", t_size, t_size);
    printf("  • 信息位列数: %d\n", base_n - t_size);
    printf("  • 码率: %.3f\n", (double)(base_n - base_m) / base_n);
    
    // ============= 配置列权重 =============
    printf("\n【第2步】 配置列权重分布（这是我们的目标！）\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    
    WD_vector col_dt;
    col_dt.n = 2;
    col_dt.wd = malloc(2 * sizeof(WD_pair));
    
    // 一半列重6
    col_dt.wd[0].weight = 6;
    col_dt.wd[0].num = base_n / 2;  // 20列
    
    // 一半列重5
    col_dt.wd[1].weight = 5;
    col_dt.wd[1].num = base_n / 2;  // 20列
    
    printf("  配置方案:\n");
    printf("    • 前%d列: 列重 = 6\n", col_dt.wd[0].num);
    printf("    • 后%d列: 列重 = 5\n", col_dt.wd[1].num);
    
    print_detail("列权重详情", &col_dt, base_n);
    
    int total_col_ones = 0;
    for (int i = 0; i < col_dt.n; i++) {
        total_col_ones += col_dt.wd[i].weight * col_dt.wd[i].num;
    }
    
    // ============= 计算行权重（自动均衡） =============
    printf("\n【第3步】 根据列权重自动计算行权重（均衡核心！）\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("  均衡原理: 矩阵中\"1\"的总数从行和列看必须相同\n");
    printf("  公式: Σ(行重 × 行数) = Σ(列重 × 列数) = %d\n", total_col_ones);
    
    WD_vector row_dt;
    row_dt.n = 2;
    row_dt.wd = malloc(2 * sizeof(WD_pair));
    
    // 计算平均行重
    int avg_row_weight = total_col_ones / base_m;
    int remainder = total_col_ones % base_m;
    
    printf("\n  计算过程:\n");
    printf("    • 总\"1\"数 = %d (来自列权重)\n", total_col_ones);
    printf("    • 总行数 = %d\n", base_m);
    printf("    • 平均行重 = %d ÷ %d = %d ... 余%d\n", 
           total_col_ones, base_m, avg_row_weight, remainder);
    printf("\n  分配策略:\n");
    printf("    • %d行: 行重 = %d (基准权重)\n", 
           base_m - remainder, avg_row_weight);
    printf("    • %d行: 行重 = %d (多分配1个以消化余数)\n", 
           remainder, avg_row_weight + 1);
    
    // 配置行重
    row_dt.wd[0].weight = avg_row_weight;
    row_dt.wd[0].num = base_m - remainder;
    
    row_dt.wd[1].weight = avg_row_weight + 1;
    row_dt.wd[1].num = remainder;
    
    print_detail("行权重详情", &row_dt, base_m);
    
    int total_row_ones = 0;
    for (int i = 0; i < row_dt.n; i++) {
        total_row_ones += row_dt.wd[i].weight * row_dt.wd[i].num;
    }
    
    // ============= 验证均衡性 =============
    printf("\n【第4步】 验证行列权重均衡性\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("  列视角统计: %d个\"1\"\n", total_col_ones);
    printf("  行视角统计: %d个\"1\"\n", total_row_ones);
    printf("  差值: %d\n", abs(total_col_ones - total_row_ones));
    
    if (total_col_ones == total_row_ones) {
        printf("\n  ✓✓✓ 均衡检查通过！行列权重完美匹配！✓✓✓\n");
    } else {
        printf("\n  ✗✗✗ 均衡检查失败！需要调整配置！✗✗✗\n");
    }
    
    // ============= 实际代码中的应用 =============
    printf("\n【第5步】 在Gen_LDPC_Codes.c中的应用示例\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("```c\n");
    printf("// 在main函数中配置\n");
    printf("WD_vector *col_dt = chk_alloc(1, sizeof(WD_vector));\n");
    printf("col_dt->n = 2;\n");
    printf("col_dt->wd = chk_alloc(2, sizeof(WD_pair));\n\n");
    
    printf("// 配置列权重：一半6，一半5\n");
    printf("col_dt->wd[0].weight = 6;\n");
    printf("col_dt->wd[0].num = %d;  // 一半列数\n", base_n/2);
    printf("col_dt->wd[1].weight = 5;\n");
    printf("col_dt->wd[1].num = %d;  // 另一半列数\n\n", base_n/2);
    
    printf("// 计算总\"1\"数\n");
    printf("int ones = 0;\n");
    printf("for (int i = 0; i < col_dt->n; i++)\n");
    printf("    ones += col_dt->wd[i].weight * col_dt->wd[i].num;\n");
    printf("// ones = %d\n\n", total_col_ones);
    
    printf("// 配置行权重：自动均衡\n");
    printf("WD_vector *row_dt = chk_alloc(1, sizeof(WD_vector));\n");
    printf("row_dt->n = 2;\n");
    printf("row_dt->wd = chk_alloc(2, sizeof(WD_pair));\n\n");
    
    printf("int m = %d;\n", base_m);
    printf("row_dt->wd[0].weight = ones / m;  // = %d\n", avg_row_weight);
    printf("row_dt->wd[0].num = m - (ones %% m);  // = %d行\n", base_m - remainder);
    printf("row_dt->wd[1].weight = ones / m + 1;  // = %d\n", avg_row_weight + 1);
    printf("row_dt->wd[1].num = ones %% m;  // = %d行\n", remainder);
    printf("```\n");
    
    // ============= 更多场景 =============
    printf("\n【额外场景】 其他列重分布示例\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    
    // 场景A: 1/3列重7, 1/3列重6, 1/3列重5
    printf("\n场景A: 三种列重均匀分布\n");
    printf("  配置: 1/3列重7 + 1/3列重6 + 1/3列重5\n");
    int ones_A = 7*(base_n/3) + 6*(base_n/3) + 5*(base_n/3);
    printf("  总\"1\"数: %d\n", ones_A);
    printf("  平均行重: %d (余%d)\n", ones_A/base_m, ones_A%base_m);
    
    // 场景B: 递减分布
    printf("\n场景B: 递减权重分布\n");
    printf("  配置: 50%%列重8 + 30%%列重6 + 20%%列重4\n");
    int ones_B = 8*(base_n/2) + 6*(base_n*3/10) + 4*(base_n/5);
    printf("  总\"1\"数: %d\n", ones_B);
    printf("  平均行重: %d (余%d)\n", ones_B/base_m, ones_B%base_m);
    
    // ============= 关键要点总结 =============
    printf("\n【关键要点总结】\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("1. 列权重可以任意指定分布（这是设计的自由度）\n");
    printf("2. 行权重必须通过总\"1\"数除以行数来自动计算\n");
    printf("3. 使用除法和取余实现行权重的尽可能均匀分布\n");
    printf("4. 核心约束: Σ(列重×列数) = Σ(行重×行数)\n");
    printf("5. 实际代码中还需考虑T矩阵的特殊处理\n");
    
    printf("\n【代码在gen_ldpc.c中的实现】\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("• 第38-56行: 检查行列\"1\"总数是否相等\n");
    printf("• 第100-106行: 反转列权重分布顺序（从高权重列开始）\n");
    printf("• 第122-145行: 初始化行计数器（row_cnt数组）\n");
    printf("• 第147-158行: 初始化列计数器（col_cnt数组）\n");
    printf("• 第180+行: 逐列填充矩阵，动态消耗row_cnt和col_cnt\n");
    
    printf("\n");
    printf("╔═══════════════════════════════════════════════════╗\n");
    printf("║              演示完成！                           ║\n");
    printf("╚═══════════════════════════════════════════════════╝\n\n");
    
    // 清理
    free(col_dt.wd);
    free(row_dt.wd);
    
    return 0;
}
