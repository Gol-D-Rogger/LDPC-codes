#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "gen_ldpc.h"
#include "alloc.h"
#include "intio.h"
#include "gen_ldpc.h"
#include "open.h"
#include "Gen_LDPC_Codes.h"

// To generate irregular LDPC according to given row and column weight distributions
int main(int argc, char **argv)
{
    read_arg(argc, argv);

    for (int file_num=1; file_num<=30; file_num++)
    {
        char filename[256];
        char cyc_filename[256];
        char mask_filename[256];
        int wt;
        int g;

        int nf = 1;

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
        sprintf(filename, "./output/%dx%d/matrix/LDPC_%dx%dex%d_w%d_dense%d_QC_H_%d_%d.txt", m-1, n-1, m-1, n-1, ex_factor, wt, g, file_num, phase);
        sprintf(mask_filename, "./output/%dx%d/mask_matrix/LDPC_%dx%dex%d_w%d_dense%d_QC_H_pad%d_mask_%d_%d.txt", m-1, n-1, m-1, n-1, ex_factor, wt, g, pad_bit, file_num, phase);
        sprintf(cyc_filename, "./output/%dx%d/cycle_record/test_%dx%d_%d_%d.txt", m-1, n-1, m-1, n-1, file_num, phase);


        // specify weight distribution

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
    
        gen_ldpc_files(filename, mask_filename, cyc_filename, nf, row_dt, col_dt, ex_factor, pad_bit, g, 4, 4, 6);
        
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

