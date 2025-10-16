#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#include "gen_ldpc.h"
#include "alloc.h"
#include "intio.h"
#include "gen_ldpc.h"
#include "open.h"
#include "Gen_LDPC_Codes.h"

// 递归创建目录（等价于 `mkdir -p`），成功或目录已存在返回 0
static int mkdir_p(const char *path)
{
    if (path == NULL || *path == '\0') return 0;

    char tmp[512];
    size_t len = strlen(path);
    if (len >= sizeof(tmp)) {
        fprintf(stderr, "mkdir_p: path too long: %s\n", path);
        return -1;
    }
    strcpy(tmp, path);

    for (char *p = tmp + 1; *p; p++)
    {
        if (*p == '/')
        {
            *p = '\0';
            if (mkdir(tmp, 0777) != 0 && errno != EEXIST)
            {
                fprintf(stderr, "mkdir_p: failed to create %s (errno=%d)\n", tmp, errno);
                return -1;
            }
            *p = '/';
        }
    }
    if (mkdir(tmp, 0777) != 0 && errno != EEXIST)
    {
        fprintf(stderr, "mkdir_p: failed to create %s (errno=%d)\n", tmp, errno);
        return -1;
    }
    return 0;
}

// To generate irregular LDPC according to given row and column weight distributions
int main(int argc, char **argv)
{
    read_arg(argc, argv);

    for (int file_num=1; file_num<=1; file_num++)
    {
        char filename[256];
        char cyc_filename[256];
        char mask_filename[256];
        int wt;
        int g;

        int nf = 1;
        int agile = 4;

        int m = base_matrix_n - base_matrix_k;
        int n = base_matrix_n;

        int ex_factor = 256;
        int pad_bit = 256;

        wt = 6;
        g = 6;
        m += 1;
        n += 1;

        int i, j, ones;
        WD_vector *row_dt, *col_dt;

        srand((unsigned)time(NULL));
        printf("################################\n");
        printf("# LDPC codes configuration: \n");
        printf("# M = %d; N= %d; P= %d; pad_bit= %d \n", m-1, n-1, ex_factor, pad_bit);
        printf("# Max. column weight = %d\n", wt);
        printf("################################\n");
        // 输出根目录支持：命令行 argv[4] 或环境变量 GENLDPC_OUT_DIR，默认 ./output
        const char *out_base_env = getenv("GENLDPC_OUT_DIR");
        const char *out_base = out_base_env ? out_base_env : (argc>=5 ? argv[4] : "./output");

        // 确保输出目录存在（等价 mkdir -p）
        char base_dir[512];
        char matrix_dir[512];
        char mask_dir[512];
        char cycle_dir[512];
        snprintf(base_dir, sizeof(base_dir), "%s/%dx%d", out_base, m-1, n-1);
        snprintf(matrix_dir, sizeof(matrix_dir), "%s/matrix", base_dir);
        snprintf(mask_dir, sizeof(mask_dir), "%s/mask_matrix", base_dir);
        snprintf(cycle_dir, sizeof(cycle_dir), "%s/cycle_record", base_dir);
        if (mkdir_p(matrix_dir) != 0 || mkdir_p(mask_dir) != 0 || mkdir_p(cycle_dir) != 0)
        {
            fprintf(stderr, "Error: cannot create output directories under %s\n", base_dir);
            exit(1);
        }

        sprintf(filename, "%s/%dx%d/matrix/LDPC_%dx%dex%d_w%d_dense%d_QC_H_%d_%d.txt",
                out_base, m-1, n-1, m-1, n-1, ex_factor, wt, g, file_num, phase);
        sprintf(mask_filename, "%s/%dx%d/mask_matrix/LDPC_%dx%dex%d_w%d_dense%d_QC_H_pad%d_mask_%d_%d.txt",
                out_base, m-1, n-1, m-1, n-1, ex_factor, wt, g, pad_bit, file_num, phase);
        sprintf(cyc_filename, "%s/%dx%d/cycle_record/test_%dx%d_%d_%d.txt",
                out_base, m-1, n-1, m-1, n-1, file_num, phase);


        // specify weight distribution
        if (agile==3)
        {
            row_dt = chk_alloc(1, sizeof(WD_vector));
            col_dt = chk_alloc(1, sizeof(WD_vector));
            col_dt->wd = chk_alloc(3, sizeof(WD_pair));
            row_dt->wd = chk_alloc(2, sizeof(WD_pair));
            col_dt->n = 3;
            row_dt->n = 2;

            ones = (n*wt - 1);
            col_dt->wd[0].weight = wt;
            col_dt->wd[0].num = base_matrix_k + 1;
            col_dt->wd[1].weight = wt - 1;
            col_dt->wd[1].num = 1;
            col_dt->wd[2].weight = wt;
            col_dt->wd[2].num = n - (base_matrix_k + 1 + 1);

            row_dt->wd[0].weight = (ones / m);
            row_dt->wd[0].num = m - (ones % m);
            row_dt->wd[1].weight = (ones / m) + 1;
            row_dt->wd[1].num = ones % m;
        }
        else if (agile==4)
        {
            row_dt = chk_alloc(1, sizeof(WD_vector));
            col_dt = chk_alloc(1, sizeof(WD_vector));
            col_dt->wd = chk_alloc(3, sizeof(WD_pair));
            row_dt->wd = chk_alloc(3, sizeof(WD_pair));
            col_dt->n = 3;
            row_dt->n = 3;

            ones = (n*wt - 1) -1 ;
            col_dt->wd[0].weight = wt;
            col_dt->wd[0].num = base_matrix_k + 1;
            col_dt->wd[1].weight = wt - 1;
            col_dt->wd[1].num = 1;
            col_dt->wd[2].weight = wt;
            col_dt->wd[2].num = n - (base_matrix_k + 1 + 1);

            row_dt->wd[0].weight = (ones / (m - 1));
            row_dt->wd[0].num = (m-1) - (ones % (m - 1));
            row_dt->wd[1].weight = (ones / (m - 1)) + 1;
            row_dt->wd[1].num = ones % (m - 1);
            row_dt->wd[2].weight = 1;
            row_dt->wd[2].num = 1;
        }
    
        gen_ldpc_files(filename, mask_filename, cyc_filename, nf, row_dt, col_dt, ex_factor, pad_bit, g, 4, 4, 6);
        // gen4_ldpc_files(filename, nf, row_dt, col_dt, ex_factor, g, 4, 4, 6);
        
        free(col_dt->wd);
        free(row_dt->wd);
        free(row_dt);
        free(col_dt);
    }
}

void read_arg(int argc, char **argv)
{
    sscanf(argv[1], "%d", &base_matrix_n);
    sscanf(argv[2], "%d", &base_matrix_k);
    sscanf(argv[3], "%d", &phase);
}
