#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
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

static void die_if_truncated(int n, size_t size, const char *label)
{
    if (n < 0 || (size_t)n >= size)
    {
        fprintf(stderr, "Error: %s path too long\n", label);
        exit(1);
    }
}

// To generate irregular LDPC according to given row and column weight distributions
int main(int argc, char **argv)
{
    read_arg(argc, argv);

    for (int file_num=1; file_num<=1; file_num++)
    {
        char filename[256];
        char occupied_filename[256];
        char fade_filename[256];
        char cyc_filename[256];
        int wt;
        int g;

        int nf = 1;
        int agile = 5;

        int m = base_matrix_n - base_matrix_k;
        int n = base_matrix_n;

        int ex_factor = 512;
        int pad_bit = 256;

        if (agile == 5)
        {
            wt = 4;
            g = 5;
        }

        int i, j, ones;
        WD_vector *row_dt, *col_dt;

        unsigned int base = (unsigned int)time(NULL);
        unsigned int pid = (unsigned int)getpid();
        const char *job_idx_env = getenv("LSB_JOBINDEX");
        unsigned int job_idx = job_idx_env ? (unsigned int)atoi(job_idx_env) : 0;
        unsigned int seed = base ^ (pid << 16) ^ job_idx ^ (unsigned int)phase;
        srand(seed);

        printf("################################\n");
        printf("# LDPC codes configuration: \n");
        printf("# M = %d; N= %d; P= %d; pad_bit= %d \n", m, n, ex_factor, pad_bit);
        printf("# Max. column weight = %d\n", wt);
        printf("################################\n");
        // 输出根目录支持：命令行 argv[4] 或环境变量 GENLDPC_OUT_DIR，默认 ./output
        const char *out_base_env = getenv("GENLDPC_OUT_DIR");
        const char *out_base = out_base_env ? out_base_env : (argc>=5 ? argv[4] : "./output");

        // 确保输出目录存在（等价 mkdir -p）
        char base_dir[512];
        char matrix_dir[512];
        char fade_dir[512];
        char occupied_dir[512];
        char cycle_dir[512];
        die_if_truncated(snprintf(base_dir, sizeof(base_dir), "%s/%dx%d", out_base, m, n),
                         sizeof(base_dir), "base_dir");
        die_if_truncated(snprintf(matrix_dir, sizeof(matrix_dir), "%s/matrix", base_dir),
                         sizeof(matrix_dir), "matrix_dir");
        die_if_truncated(snprintf(fade_dir, sizeof(fade_dir), "%s/fade_matrix", base_dir),
                         sizeof(fade_dir), "fade_dir");
        die_if_truncated(snprintf(occupied_dir, sizeof(occupied_dir), "%s/occupied_matrix", base_dir),
                         sizeof(occupied_dir), "occupied_dir");
        die_if_truncated(snprintf(cycle_dir, sizeof(cycle_dir), "%s/cycle_record", base_dir),
                         sizeof(cycle_dir), "cycle_dir");
        if (mkdir_p(matrix_dir) != 0 || mkdir_p(fade_dir) != 0 || mkdir_p(occupied_dir) != 0 || mkdir_p(cycle_dir) != 0)
        {
            fprintf(stderr, "Error: cannot create output directories under %s\n", base_dir);
            exit(1);
        }

        sprintf(filename, "%s/%dx%d/matrix/LDPC_%dx%dex%d_w%d_dense%d_QC_H_%d_%d.txt",
                out_base, m, n, m, n, ex_factor, wt, g, file_num, phase);
        sprintf(fade_filename, "%s/%dx%d/fade_matrix/LDPC_%dx%dex%d_w%d_dense%d_fade_%d_%d.txt",
                out_base, m, n, m, n, ex_factor, wt, g, file_num, phase);
        sprintf(occupied_filename, "%s/%dx%d/occupied_matrix/LDPC_%dx%dex%d_w%d_dense%d_occupied_%d_%d.txt",
                out_base, m, n, m, n, ex_factor, wt, g, file_num, phase);
        sprintf(cyc_filename, "%s/%dx%d/cycle_record/LDPC_%dx%dex%d_w%d_dense%d_cycle_%d_%d.txt",
                out_base, m, n, m, n, ex_factor, wt, g, file_num, phase);


        // specify weight distribution
        if (agile==5)
        {
            col_dt = chk_alloc(1, sizeof(WD_vector));
            row_dt = chk_alloc(1, sizeof(WD_vector));
            col_dt->wd = chk_alloc(2, sizeof(WD_pair));
            row_dt->wd = chk_alloc(3, sizeof(WD_pair));
            col_dt->n = 2;
            row_dt->n = 3;

            ones = (n*wt - 1);
            col_dt->wd[0].weight = wt;
            col_dt->wd[0].num = base_matrix_n - 1;
            col_dt->wd[1].weight = wt - 1;
            col_dt->wd[1].num = 1;

            if (m > 5)
            {
                row_dt->wd[0].weight = (ones / (m-1) + 1);
                row_dt->wd[0].num = ones%(m-1);
                row_dt->wd[1].weight =ones/(m-1);
                row_dt->wd[1].num = (m-1)-ones%(m-1);
                row_dt->wd[2].weight = 0;
                row_dt->wd[2].num = 1;
            }
            else
            {
                row_dt->wd[0].weight = ones / m + 1;
                row_dt->wd[0].num =ones % m;
                row_dt->wd[1].weight = ones / m;
                row_dt->wd[1].num = m - ones % m;
                row_dt->wd[2].weight = 0;
                row_dt->wd[2].num = 0;
            }
        }
        
        gen_ldpc_files(filename, occupied_filename, fade_filename, cyc_filename,
                       nf, row_dt, col_dt, ex_factor, pad_bit, g, 4, 4, 6);
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
