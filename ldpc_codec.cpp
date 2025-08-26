#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rand.h"
#include "mod2sparse.h"
#include "mod2dense.h"
#include "mod2convert.h"
#include "vec_op.h"
#include "finite_lib.h"
#include "transceiver.h"
#include "ldpc_codec.h"

void ldpc_packet::ldpc_rd_phck(char *pchk_file, char *drop_file)
{
    mod2entry *e;
    FILE *fp;
    FILE *fp_drop;
    int col_shift;
    qc_bm = mod2sparse_allocate(bm_m, bm_n);
    qc_hm = mod2sparse_allocate(hm_m, hm_n);

    printf("[LDPC] Porting H matrix from %s ...\n", pchk_file);

    fp = fopen(pchk_file, "r");
    fp_drop = fopen(drop_file, "r");

    if (fp == NULL || fp_drop == NULL)
    {
        printf("[LDPC] Error opening file %s or %s\n", pchk_file, drop_file);
        exit(1);
    }


    drop_col = new int *[bm_m];
    for(int i=0;i<bm_m-1;i++)
        drop_col[i] = new int [8];
    drop_col_bit_map = new int*[bm_m];
    drop_col[bm_m-1] = new int [80];
    for(int i=0;i<bm_m;i++)
    {
        drop_col_bit_map[i] = new int[bm_m];
        for(int j=0;j<bm_n;j++)
            drop_col_bit_map[i][j] = 0;
    }


    for(int i=0;i<80;i++)
        drop_col[bm_m-1][i]=-1;
    wit_drop = new int [bm_m];
    for(int i=0;i<bm_m;i++)
        wit_drop[i]=0;
    for(int j=0;j<8;j++)
    {
        for(int i=0;i<bm_m-1;i++)
        {
            fscanf(fp_drop,"%d",&drop_col[i][j]);
            if(drop_col[i][j] != -1)
            {
                drop_col_bit_map[i][drop_col[i][j]] = 1;
                drop_col_bit_map[bm_m-1][drop_col[i][j]] = 1;
                wit_drop[i]++;
                drop_col[bm_m-1][wit_drop[bm_m-1]++]=drop_col[i][j];
            }
        }
    }
    fclose(fp_drop);


    int remain_length = cir_sz-pad_bit;
    for(int i=0;i<bm_m-1;i++)
    {
        for(int j=0;j<bm_n;j++)
        {
            fscanf(fp,"%d",&col_shift);
            if((i<tm_sz) && (j>(bm_n - tm_sz-2)) && j!=bm_n-1)
            {
                if(((i!=(j+tm_sz-bm_n+1))&&(col_shift>=0)) || ((i=(j+tm_sz-bm_n+1)) &&(col_shift!=0)))
                {
                    //submatrix T is not identity matrix
                    exit(0);
                }
            }


            if(col_shift>=0)
            {
                if(pad_bit == cir_sz && drop_col_bit_map[i][j])
                {}
                else
                {
                    e = mod2sparse_insert(qc_bm,i,j);
                    e->shift = col_shift;
                }
                if(drop_col_bit_map[i][j])
                {
                    int rest_start = ((j-i)/(bm_m-1)+1)*pad_bit %cir_sz;
                    if(remain_length+rest_start<=cir_sz)
                    {
                      for (int k = rest_start; k < rest_start + remain_length; k++)
                        mod2sparse_insert(qc_hm, i * cir_sz + k,
                                          j * cir_sz +
                                              (k + col_shift) % cir_sz);
                    }
                    else
                    {
                        for(int k =rest_start;k<cir_sz;k++)
                            mod2sparse_insert(qc_hm,i*cir_sz+k,j*cir_sz + (k+col_shift)%cir_sz);
                        for(int k =0 ;k<remain_length-cir_sz+rest_start;k++)
                            mod2sparse_insert(qc_hm,i*cir_sz+k,j*cir_sz + (k+col_shift)%cir_sz);
                    }


                    if(pad_bit !=0)
                    {
                        int mask_row_ind = (j-i)/(bm_m-1);
                        e = mod2sparse_insert(qc_bm,(bm_m-1),j);
                        col_shift = (col_shift + mask_row_ind * pad_bit) % cir_sz;
                        e->shift = col_shift;
                        for(int k=0; k< pad_bit;k++)
                            mod2sparse_insert(qc_hm, (bm_m-1) * cir_sz + k, j*cir_sz + (k+col_shift) %cir_sz);
                    }
                }


                else
                {
                    for(int k = 0;k<cir_sz;k++)
                        mod2sparse_insert(qc_bm, i*cir_sz+k,j*cir_sz+(k+col_shift) % cir_sz);
                }
            }
        }
    }


    if(pad_bit != 0)
    {
        e = mod2sparse_insert(qc_bm,(bm_m-1),bm_n-1);
        e->shift = 0;
        for(int k=0; k< cir_sz;k++)
            mod2sparse_insert(qc_hm, (bm_m-1) * cir_sz + k, (bm_n-1)*cir_sz+k);
    }


    fclose(fp);
    printf("[LDPC] H matrix porting ready!\n");
}


// Generate G matrices from H
void ldpc_packet::ldpc_gen_gm()
{
    FILE *fp;
    mod2sparse *qc_ac, *qc_bd, *qc_te, *qc_fall;
    mod2sparse *qc_exb, *qc_f;
    mod2dense *qc_f_d, *qc_fi_d;
    int *qc_a_cols, *qc_b_cols, *qc_t_cols, *qc_fall_cols;
    int *qc_a_rows, *qc_c_rows, *qc_g_rows;      

    qc_ac = mod2sparse_allocate(hm_m, hm_k);
    qc_bd = mod2sparse_allocate(hm_m, (bm_m-tm_sz-1)*cir_sz);
    qc_te = mod2sparse_allocate(hm_m, tm_sz*cir_sz);
    qc_fall = mod2sparse_allocate(hm_m, cir_sz);

    qc_a_cols = (int *)calloc((bm_n-bm_m)*cir_sz, sizeof(*qc_a_cols));
    qc_b_cols = (int *)calloc((bm_m-tm_sz-1)*cir_sz, sizeof(*qc_b_cols));
    qc_t_cols = (int *)calloc((tm_sz*cir_sz), sizeof(*qc_t_cols));
    qc_fall_cols = (int *)calloc(cir_sz, sizeof(*qc_fall_cols));
    qc_a_rows = (int *)calloc(tm_sz*cir_sz, sizeof(*qc_a_rows));
    qc_c_rows = (int *)calloc((bm_m-tm_sz-1)*cir_sz, sizeof(*qc_c_rows));
    qc_g_rows = (int *)calloc(cir_sz, sizeof(*qc_g_rows));

    printf("[LDPC] Generating A/B/C/D/E matrices from H ...\n");

    for (int i=0; i<bm_k*cir_sz; i++)
        qc_a_cols[i] = i;
    for (int i=0; i<(bm_m-tm_sz-1)*cir_sz; i++)
        qc_b_cols[i] = i + bm_k*cir_sz;
    for (int i=0; i<(tm_sz*cir_sz); i++)
        qc_t_cols[i] = i + (bm_m-tm_sz-1)*cir_sz;
    for (int i=0; i<cir_sz; i++)
        qc_fall_cols[i] = i + (bm_n-1)*cir_sz;

    for (int i=0; i<tm_sz*cir_sz; i++)
        qc_a_rows[i] = i;
    for (int i=0; i<(bm_m-tm_sz-1)*cir_sz; i++)
        qc_c_rows[i] = i + tm_sz*cir_sz;
    for (int i=0; i<cir_sz; i++)
        qc_g_rows[i] = i + (bm_m-1)*cir_sz;

    // split column first
    mod2sparse_copycols(qc_hm, qc_ac, qc_a_cols);
    mod2sparse_copycols(qc_hm, qc_bd, qc_b_cols);
    mod2sparse_copycols(qc_hm, qc_te, qc_t_cols);
    mod2sparse_copycols(qc_hm, qc_fall, qc_fall_cols);

    // then split row
    mod2sparse_copyrows(qc_ac, qc_a, qc_a_rows);
    mod2sparse_copyrows(qc_ac, qc_c, qc_c_rows);
    mod2sparse_copyrows(qc_ac, qc_g, qc_g_rows);
    mod2sparse_copyrows(qc_bd, qc_b, qc_a_rows);
    mod2sparse_copyrows(qc_bd, qc_d, qc_c_rows);
    mod2sparse_copyrows(qc_te, qc_e, qc_c_rows);
    mod2sparse_copyrows(qc_fall, qc_f1, qc_a_rows);
    mod2sparse_copyrows(qc_fall, qc_f2, qc_c_rows);
    printf("[LDPC] A/B/C/D/E matrices generated!\n");

#ifdef _LDPC_DEBUG_DUMP
    char MH[50] = "H_matrix.txt";
    char MA[50] = "A_matrix.txt";
    char MB[50] = "B_matrix.txt";
    char MC[50] = "C_matrix.txt";
    char MD[50] = "D_matrix.txt";
    char ME[50] = "E_matrix.txt";

    printf("[LDPC] Dump H matrix to file %s\n", MH);
    fp = fopen(MH, "w");
    mod2sprase_print(fp, qc,hm);
    fclose(fp);
    printf("[LDPC] Dump A matrix to file %s\n", MA);
    fp = fopen(MA, "w");
    mod2sprase_print(fp, qc_a);
    fclose(fp);
    printf("[LDPC] Dump B matrix to file %s\n", MB);
    fp = fopen(MB, "w");
    mod2sprase_print(fp, qc_b);
    fclose(fp);
    printf("[LDPC] Dump C matrix to file %s\n", MC);
    fp = fopen(MC, "w");
    mod2sprase_print(fp, qc_c);
    fclose(fp);
    printf("[LDPC] Dump D matrix to file %s\n", MD);
    fp = fopen(MD, "w");
    mod2sprase_print(fp, qc_d);
    fclose(fp);
    printf("[LDPC] Dump E matrix to file %s\n", ME);
    fp = fopen(ME, "w");
    mod2sprase_print(fp, qc_e);
    fclose(fp);
#endif    

    // generate inverse F matrix (F=E*B+D)
    printf("[LDPC] Generating inverse F matrix from H matrix ...\n");

    qc_exb = mod2sparse_allocate((bm_m-tm_sz-1)*cir_sz, (bm_m-tm_sz-1)*cir_sz);
    qc_f = mod2sparse_allocate((bm_m-tm_sz-1)*cir_sz, (bm_m-tm_sz-1)*cir_sz);

    mod2sparse_multiply(qc_e, qc_b, qc_exb);
    mod2sparse_add(qc_exb, qc_d, qc_f);

    qc_f_d = mod2dense_allocate((bm_m-tm_sz-1)*cir_sz, (bm_m-tm_sz-1)*cir_sz);
    qc_fi_d = mod2dense_allocate((bm_m-tm_sz-1)*cir_sz, (bm_m-tm_sz-1)*cir_sz);

    mod2sparse_to_dense(qc_f, qc_f_d);
    mod2dense_invert(qc_f_d, qc_fi_d);
    mod2dense_to_sparse(qc_fi_d, qc_f);

    printf("[LDPC] Inverse F matrix generated!\n");

#ifdef _LDPC_DEBUG_DUMP
    char MFi[50] ;
    sprintf(MFi, "./output/LDPC_%dx%dwx%d_w%d_dense6_Fi_matrix.txt", bm_m, bm_n, cir_sz, col_wt);
    printf("[LDPC] Dump inverse F matrix to file %s\n", MFi);
    fp = fopen(MFi, "w");
    mod2sprase_print(fp, qc_fi);
    fclose(fp);
#endif

    // free up all internals
    mod2sparse_free(qc_ac);
    mod2sparse_free(qc_bd);
    mod2sparse_free(qc_te);
    mod2sparse_free(qc_fall);
    mod2sparse_free(qc_exb);
    mod2sparse_free(qc_f);
    mod2dense_free(qc_f_d);
    mod2dense_free(qc_fi_d);
    free(qc_a_cols);
    free(qc_b_cols);
    free(qc_t_cols);
    free(qc_fall_cols);
    free(qc_a_rows);
    free(qc_c_rows);
    free(qc_g_rows);
} // ldpc_gen_gm


// LDPC code configuration
void ldpc_packet::ldpc_config(int m, int n, int sc, int pdbit, int st, int wt, int always_on, char *pchk_file, char *drop_file)
{
    // QC matrix sizes
    bm_m = m;
    bm_n = n;
    bm_k = n - m;
    tm_sz = st;
    cir_sz = sc;
    awon = always_on;

    col_wt = wt;
    pad_bit = pdbit;
    hm_m = m * sc;
    hm_n = n * sc;
    hm_k = hm_n - hm_m;
    pad_len = hm_k - info_len;

    printf("[LDPC] Configuring LDPC code with m=%d, n=%d, sc=%d, st=%d, wt=%d\n", 
           bm_m, bm_n, cir_sz, tm_sz, col_wt);

    // read parity check matrix
    ldpc_rd_phck(pchk_file, drop_file);

    // G matrices
    qc_a = mod2sparse_allocate(tm_sz * cir_sz, (bm_k) * cir_sz);
    qc_b = mod2sparse_allocate(tm_sz * cir_sz, (bm_m - tm_sz - 1) * cir_sz);
    qc_c = mod2sparse_allocate((bm_m - tm_sz - 1) * cir_sz, bm_k * cir_sz);
    qc_d = mod2sparse_allocate((bm_m - tm_sz - 1) * cir_sz, (bm_m - tm_sz - 1) * cir_sz);
    qc_e = mod2sparse_allocate((bm_m - tm_sz - 1) * cir_sz, (tm_sz) * cir_sz);
    qc_g = mod2sparse_allocate(cir_sz, bm_k * cir_sz);
    qc_f1 = mod2sparse_allocate(tm_sz * cir_sz, cir_sz);
    qc_f2 = mod2sparse_allocate((bm_m - tm_sz - 1) * cir_sz, cir_sz);
    qc_fi = mod2sparse_allocate((bm_m - tm_sz - 1) * cir_sz, (bm_m - tm_sz - 1) * cir_sz);

    // generate G matrices
    ldpc_gen_gm();

#ifdef _LDPC_DEBUG_DUMP
        rdec_cmem_cont_thrshd = 10;
        rdec_hdmem_cont_thrshd = 6;

        // print_hm();
#endif
} // ldpc_config


// LDPC decoder config
void ldpc_packet::ldpc_dec_config(int max_fdec_itr, int max_tbfdec_itr, int fdec_col_skip, int max_ldec_itr, float dec_alpha, int fin_mode, int fin_q_num, int fin_r_num, int fin_f_num, int fp_flg)
{
    // syndrome weight
    init_synd_wt_min = hm_m;
    init_synd_wt_max = 0;

    reg_fp_flg = fp_flg;

    // layer config
    if (max_ldec_itr>0)
    {
        ldec_max_itr = max_ldec_itr;
        ldec_early_term_en = 1;
    }
    else
    {
        ldec_max_itr = -1*max_ldec_itr;
        ldec_early_term_en = 0;
    }

    alpha = dec_alpha;

    //initialize quantization
    finite_mode = fin_mode;
    finite_q_num = fin_q_num;
    finite_r_num = fin_r_num;
    finite_c_num = finite_r_num;
    finite_f_num = fin_f_num;

    int tmp = pow(2, finite_f_num);
    finite_q_max = (pow(2, finite_q_num-1) -1)/tmp;
    finite_q_min = -(pow(2, finite_q_num-1) -1)/tmp;
    finite_r_max = (pow(2, finite_r_num-1) -1)/tmp;
    finite_r_min = -(pow(2, finite_r_num-1) -1)/tmp;
    finite_c_max = (pow(2, finite_c_num-1) -1)/tmp;
    finite_c_min = 0;

#ifdef _LDPC_DEBUG_DUMP
    printf("[LDPC DEBUG] Q MSG: %d/%d, max %f, min %f \n", finite_q_num, finite_q_max, finite_q_min);
    printf("[LDPC DEBUG] R MSG: %d/%d, max %f, min %f \n", finite_r_num, finite_r_max, finite_r_min);
    printf("[LDPC DEBUG] C MSG: %d/%d, max %f, min %f \n", finite_c_num-1, finite_c_max, finite_c_min);
#endif    

    // BF config
    if (max_fdec_itr>0)
    {
        fdec_max_itr = max_fdec_itr;
        fdec_early_term_en = 1;
    }
    else
    {
        fdec_max_itr = -1*max_fdec_itr;
        fdec_early_term_en = 0;
    }

    if (max_tbfdec_itr>0)
    {
        tbfdec_max_itr = max_tbfdec_itr;
        fdec_early_term_en = 1;
    }
    else
    {
        tbfdec_max_itr = -1*max_tbfdec_itr;
        fdec_early_term_en = 0;
    }

    col_skip_itr = fdec_col_skip;

    flp_thrshd0 = (int *)calloc(fdec_max_itr, sizeof(*flp_thrshd0));
    flp_thrshd1 = (int *)calloc(fdec_max_itr, sizeof(*flp_thrshd1));

    flp_thrshd0[0] = 6;
    flp_thrshd0[1] = 5;
    flp_thrshd0[2] = 5;
    flp_thrshd0[3] = 5;
    flp_thrshd0[4] = 5;
    flp_thrshd0[5] = 5;
    flp_thrshd0[6] = 5;
    flp_thrshd0[7] = 5;
    for (int i=8; i<fdec_max_itr;i++)
        flp_thrshd0[i] = 4;

    flp_thrshd1[0] = 3;
    flp_thrshd1[1] = 1;
    flp_thrshd1[2] = 1;
    flp_thrshd1[3] = 1;
    flp_thrshd1[4] = 1;
    flp_thrshd1[5] = 1;
    flp_thrshd1[6] = 2;
    flp_thrshd1[7] = 4;
    for (int i=8; i<fdec_max_itr;i++)
        flp_thrshd1[i] = 2;

    tbbf_thrshd0 = (int *)calloc(128, sizeof(*tbbf_thrshd0));
    tbbf_thrshd1 = (int *)calloc(128, sizeof(*tbbf_thrshd1));
    int mxtbbf_thrshd0[128] = {0};
    int mxtbbf_thrshd1[128] = {0};

    mxtbbf_thrshd0[0] = 3;  mxtbbf_thrshd1[0] = -6;
    mxtbbf_thrshd0[1] = 2;  mxtbbf_thrshd1[1] = 1;
    mxtbbf_thrshd0[2] = 2;  mxtbbf_thrshd1[2] = -1;
    mxtbbf_thrshd0[3] = 2;  mxtbbf_thrshd1[3] = -3;
    mxtbbf_thrshd0[4] = 2;  mxtbbf_thrshd1[4] = 3;
    mxtbbf_thrshd0[5] = 2;  mxtbbf_thrshd1[5] = -3;
    mxtbbf_thrshd0[6] = 3;  mxtbbf_thrshd1[6] = 3;
    mxtbbf_thrshd0[7] = 2;  mxtbbf_thrshd1[7] = -1;
    mxtbbf_thrshd0[8] = 1;  mxtbbf_thrshd1[8] = -3;
    mxtbbf_thrshd0[9] = 2;  mxtbbf_thrshd1[9] = 3;
    mxtbbf_thrshd0[10] = 2;  mxtbbf_thrshd1[10] = -3;
    mxtbbf_thrshd0[11] = 2;  mxtbbf_thrshd1[11] = 5;
    mxtbbf_thrshd0[12] = 1;  mxtbbf_thrshd1[12] = -1;
    mxtbbf_thrshd0[13] = 1;  mxtbbf_thrshd1[13] = -3;
    mxtbbf_thrshd0[14] = 2;  mxtbbf_thrshd1[14] = 3;
    mxtbbf_thrshd0[15] = 1;  mxtbbf_thrshd1[15] = -1;
    mxtbbf_thrshd0[16] = 1;  mxtbbf_thrshd1[16] = 0;
    mxtbbf_thrshd0[17] = 1;  mxtbbf_thrshd1[17] = -3;
    mxtbbf_thrshd0[18] = 1;  mxtbbf_thrshd1[18] = 1;
    mxtbbf_thrshd0[19] = 1;  mxtbbf_thrshd1[19] = -1;
    mxtbbf_thrshd0[20] = 1;  mxtbbf_thrshd1[20] = 1;
    mxtbbf_thrshd0[21] = 1;  mxtbbf_thrshd1[21] = 0;
    mxtbbf_thrshd0[22] = 1;  mxtbbf_thrshd1[22] = 0;
    mxtbbf_thrshd0[23] = 1;  mxtbbf_thrshd1[23] = 1;
    mxtbbf_thrshd0[24] = 0;  mxtbbf_thrshd1[24] = 0;
    mxtbbf_thrshd0[25] = 1;  mxtbbf_thrshd1[25] = 0;
    mxtbbf_thrshd0[26] = 0;  mxtbbf_thrshd1[26] = 0;
    mxtbbf_thrshd0[27] = 0;  mxtbbf_thrshd1[27] = 0;
    mxtbbf_thrshd0[28] = 1;  mxtbbf_thrshd1[28] = 0;
    mxtbbf_thrshd0[29] = 1;  mxtbbf_thrshd1[29] = 0;
    mxtbbf_thrshd0[30] = 1;  mxtbbf_thrshd1[30] = 1;
    mxtbbf_thrshd0[31] = 1;  mxtbbf_thrshd1[31] = 1;   
    mxtbbf_thrshd0[32] = 1;  mxtbbf_thrshd1[32] = 1;
    mxtbbf_thrshd0[33] = 1;  mxtbbf_thrshd1[33] = 1;
    mxtbbf_thrshd0[34] = 1;  mxtbbf_thrshd1[34] = 1;
    mxtbbf_thrshd0[35] = 1;  mxtbbf_thrshd1[35] = 1;
    mxtbbf_thrshd0[36] = 1;  mxtbbf_thrshd1[36] = 1;
    for (int i=0; i<37; i++)
    {
        tbbf_thrshd0[i] = mxtbbf_thrshd0[i];
        tbbf_thrshd1[i] = mxtbbf_thrshd1[i];
    }
    for (int i=37; i<128; i++)
    {
        tbbf_thrshd0[i] = 1;
        tbbf_thrshd1[i] = 1;
    }

    // special iteration 
    // break some deadlocks that would degrade performance
    // especially for matrix 143x14wt6

    flp_thrshd0_s = (int*)calloc(fdec_max_itr, sizeof(*flp_thrshd0_s));
    flp_thrshd0_w = (int*)calloc(fdec_max_itr, sizeof(*flp_thrshd0_w));
    flp_thrshd1_s = (int*)calloc(fdec_max_itr, sizeof(*flp_thrshd1_s));
    flp_thrshd1_w = (int*)calloc(fdec_max_itr, sizeof(*flp_thrshd1_w));

    flp_thrshd0_s[0] = 6;
    flp_thrshd0_s[1] = 6;
    flp_thrshd0_s[2] = 6;
    flp_thrshd0_s[3] = 6;
    flp_thrshd0_s[4] = 6;
    flp_thrshd0_s[5] = 6;
    flp_thrshd0_s[6] = 6;
    flp_thrshd0_s[7] = 6;
    for (int i=8; i<fdec_max_itr;i++)
        flp_thrshd0_s[i] = 5;

    flp_thrshd0_w[0] = 6;
    flp_thrshd0_w[1] = 5;
    flp_thrshd0_w[2] = 5;
    flp_thrshd0_w[3] = 5;
    flp_thrshd0_w[4] = 5;
    flp_thrshd0_w[5] = 5;
    flp_thrshd0_w[6] = 5;
    flp_thrshd0_w[7] = 5;
    for (int i=8; i<fdec_max_itr;i++)
        flp_thrshd0[i] = 4;

    flp_thrshd1_s[0] = 3;
    flp_thrshd1_s[1] = 2;
    flp_thrshd1_s[2] = 2;
    flp_thrshd1_s[3] = 2;
    flp_thrshd1_s[4] = 2;
    flp_thrshd1_s[5] = 2;
    flp_thrshd1_s[6] = 2;
    flp_thrshd1_s[7] = 4;
    for (int i=8; i<fdec_max_itr;i++)
        flp_thrshd1_s[i] = 3;

    flp_thrshd1_w[0] = 3;
    flp_thrshd1_w[1] = 1;
    flp_thrshd1_w[2] = 1;
    flp_thrshd1_w[3] = 1;
    flp_thrshd1_w[4] = 1;
    flp_thrshd1_w[5] = 1;
    flp_thrshd1_w[6] = 2;
    flp_thrshd1_w[7] = 4;
    for (int i=8; i<fdec_max_itr;i++)
        flp_thrshd1_w[i] = 2;

    sb_thrshd0_s0 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd0_s0));
    sb_thrshd0_s1 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd0_s1));
    sb_thrshd0_w0 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd0_w0));
    sb_thrshd0_w1 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd0_w1));
    sb_thrshd1_s0 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd1_s0));
    sb_thrshd1_s1 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd1_s1));
    sb_thrshd1_w0 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd1_w0));
    sb_thrshd1_w1 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd1_w1));

    sb_thrshd0_s0[0] = 6;
    sb_thrshd0_s0[1] = 4;
    sb_thrshd0_s0[2] = 4;
    sb_thrshd0_s0[3] = 4;
    sb_thrshd0_s0[4] = 4;
    sb_thrshd0_s0[5] = 4;
    sb_thrshd0_s0[6] = 4;
    sb_thrshd0_s0[7] = 4;
    for (int i=8; i<fdec_max_itr;i++)
        sb_thrshd0_s0[i] = 3;

    sb_thrshd0_s1[0] = 7;
    sb_thrshd0_s1[1] = 7;
    sb_thrshd0_s1[2] = 7;
    sb_thrshd0_s1[3] = 7;
    sb_thrshd0_s1[4] = 7;
    sb_thrshd0_s1[5] = 7;
    sb_thrshd0_s1[6] = 7;
    sb_thrshd0_s1[7] = 7;
    for (int i=8; i<fdec_max_itr;i++)
        sb_thrshd0_s1[i] = 7;

    sb_thrshd0_w0[0] = 3;
    sb_thrshd0_w0[1] = 3;
    sb_thrshd0_w0[2] = 3;
    sb_thrshd0_w0[3] = 3;
    sb_thrshd0_w0[4] = 3;
    sb_thrshd0_w0[5] = 3;
    sb_thrshd0_w0[6] = 3;
    sb_thrshd0_w0[7] = 3;
    for (int i=8; i<fdec_max_itr;i++)
        sb_thrshd0_s1[i] = 2;

    sb_thrshd0_w1[0] = 6;
    sb_thrshd0_w1[1] = 6;
    sb_thrshd0_w1[2] = 6;
    sb_thrshd0_w1[3] = 6;
    sb_thrshd0_w1[4] = 6;
    sb_thrshd0_w1[5] = 6;
    sb_thrshd0_w1[6] = 6;
    sb_thrshd0_w1[7] = 6;
    for (int i=8; i<fdec_max_itr;i++)
        sb_thrshd0_w1[i] = 6;

    sb_thrshd1_s0[0] = 0;
    sb_thrshd1_s0[1] = 0;
    sb_thrshd1_s0[2] = 0;
    sb_thrshd1_s0[3] = 0;
    sb_thrshd1_s0[4] = 0;
    sb_thrshd1_s0[5] = 0;
    sb_thrshd1_s0[6] = 0;
    sb_thrshd1_s0[7] = 0;
    for (int i=8; i<fdec_max_itr;i++)
        sb_thrshd1_s0[i] = 0;

    sb_thrshd1_s1[0] = 6;
    sb_thrshd1_s1[1] = 6;
    sb_thrshd1_s1[2] = 6;
    sb_thrshd1_s1[3] = 6;
    sb_thrshd1_s1[4] = 6;
    sb_thrshd1_s1[5] = 6;
    sb_thrshd1_s1[6] = 6;
    sb_thrshd1_s1[7] = 6;
    for (int i=8; i<fdec_max_itr;i++)
        sb_thrshd1_s1[i] = 6;

    sb_thrshd1_w0[0] = 0;
    sb_thrshd1_w0[1] = 0;
    sb_thrshd1_w0[2] = 0;
    sb_thrshd1_w0[3] = 0;
    sb_thrshd1_w0[4] = 0;
    sb_thrshd1_w0[5] = 0;
    sb_thrshd1_w0[6] = 0;
    sb_thrshd1_w0[7] = 0;
    for (int i=8; i<fdec_max_itr;i++)
        sb_thrshd1_w0[i] = 0;

    sb_thrshd1_w1[0] = 5;
    sb_thrshd1_w1[1] = 5;
    sb_thrshd1_w1[2] = 5;
    sb_thrshd1_w1[3] = 5;
    sb_thrshd1_w1[4] = 5;
    sb_thrshd1_w1[5] = 5;
    sb_thrshd1_w1[6] = 5;
    sb_thrshd1_w1[7] = 5;
    for (int i=8; i<fdec_max_itr;i++)
        sb_thrshd1_w1[i] = 5;

    fpd_flp_thrshd0_s = (int *)calloc(fdec_max_itr, sizeof(*fpd_flp_thrshd0_s));
    fpd_flp_thrshd0_w = (int *)calloc(fdec_max_itr, sizeof(*fpd_flp_thrshd0_w));
    fpd_flp_thrshd1_s = (int *)calloc(fdec_max_itr, sizeof(*fpd_flp_thrshd1_s));
    fpd_flp_thrshd1_w = (int *)calloc(fdec_max_itr, sizeof(*fpd_flp_thrshd1_w));

    fpd_flp_thrshd0_s[0] = 6;
    fpd_flp_thrshd0_s[1] = 6;
    fpd_flp_thrshd0_s[2] = 6;
    fpd_flp_thrshd0_s[3] = 6;
    fpd_flp_thrshd0_s[4] = 6;
    fpd_flp_thrshd0_s[5] = 6;
    fpd_flp_thrshd0_s[6] = 6;
    fpd_flp_thrshd0_s[7] = 6;
    for (int i=8; i<fdec_max_itr;i++)
        fpd_flp_thrshd0_s[i] = 5;

    fpd_flp_thrshd0_w[0] = 6;
    fpd_flp_thrshd0_w[1] = 5;
    fpd_flp_thrshd0_w[2] = 5;
    fpd_flp_thrshd0_w[3] = 5;
    fpd_flp_thrshd0_w[4] = 5;
    fpd_flp_thrshd0_w[5] = 5;
    fpd_flp_thrshd0_w[6] = 5;
    fpd_flp_thrshd0_w[7] = 5;
    fpd_flp_thrshd0_w[8] = 5;
    for (int i=9; i<fdec_max_itr;i++)
        fpd_flp_thrshd0_w[i] = 4;

    fpd_flp_thrshd1_s[0] = 3;
    fpd_flp_thrshd1_s[1] = 2;
    fpd_flp_thrshd1_s[2] = 2;
    fpd_flp_thrshd1_s[3] = 2;
    fpd_flp_thrshd1_s[4] = 2;
    fpd_flp_thrshd1_s[5] = 2;
    fpd_flp_thrshd1_s[6] = 2;
    fpd_flp_thrshd1_s[7] = 4;
    for (int i=8; i<fdec_max_itr;i++)
        fpd_flp_thrshd1_s[i] = 3;

    fpd_flp_thrshd1_w[0] = 3;
    fpd_flp_thrshd1_w[1] = 1;
    fpd_flp_thrshd1_w[2] = 1;
    fpd_flp_thrshd1_w[3] = 1;
    fpd_flp_thrshd1_w[4] = 1;
    fpd_flp_thrshd1_w[5] = 1;
    fpd_flp_thrshd1_w[6] = 2;
    fpd_flp_thrshd1_w[7] = 4;
    for (int i=8; i<fdec_max_itr;i++)
        fpd_flp_thrshd1_w[i] = 2;

    fpd_sb_thrshd0_s0 = (int *)calloc(fdec_max_itr, sizeof(*fpd_sb_thrshd0_s0));
    fpd_sb_thrshd0_s1 = (int *)calloc(fdec_max_itr, sizeof(*fpd_sb_thrshd0_s1));
    fpd_sb_thrshd0_w0 = (int *)calloc(fdec_max_itr, sizeof(*fpd_sb_thrshd0_w0));
    fpd_sb_thrshd0_w1 = (int *)calloc(fdec_max_itr, sizeof(*fpd_sb_thrshd0_w1));
    fpd_sb_thrshd1_s0 = (int *)calloc(fdec_max_itr, sizeof(*fpd_sb_thrshd1_s0));
    fpd_sb_thrshd1_s1 = (int *)calloc(fdec_max_itr, sizeof(*fpd_sb_thrshd1_s1));
    fpd_sb_thrshd1_w0 = (int *)calloc(fdec_max_itr, sizeof(*fpd_sb_thrshd1_w0));
    fpd_sb_thrshd1_w1 = (int *)calloc(fdec_max_itr, sizeof(*fpd_sb_thrshd1_w1));

    fpd_sb_thrshd0_s0[0] = 6;
    fpd_sb_thrshd0_s0[1] = 4;
    fpd_sb_thrshd0_s0[2] = 5;
    fpd_sb_thrshd0_s0[3] = 4;
    fpd_sb_thrshd0_s0[4] = 4;
    fpd_sb_thrshd0_s0[5] = 4;
    fpd_sb_thrshd0_s0[6] = 4;
    fpd_sb_thrshd0_s0[7] = 4;
    for (int i=8; i<fdec_max_itr;i++)
        fpd_sb_thrshd0_s0[i] = 4;

    fpd_sb_thrshd0_s1[0] = 7;
    fpd_sb_thrshd0_s1[1] = 7;
    fpd_sb_thrshd0_s1[2] = 7;
    fpd_sb_thrshd0_s1[3] = 7;
    fpd_sb_thrshd0_s1[4] = 7;
    fpd_sb_thrshd0_s1[5] = 7;
    fpd_sb_thrshd0_s1[6] = 7;
    fpd_sb_thrshd0_s1[7] = 7;
    fpd_sb_thrshd0_s1[8] = 7;
    fpd_sb_thrshd0_s1[9] = 7;
    for (int i=10; i<fdec_max_itr;i++)
        fpd_sb_thrshd0_s1[i] = 7;

    fpd_sb_thrshd0_w0[0] = 3;
    fpd_sb_thrshd0_w0[1] = 3;
    fpd_sb_thrshd0_w0[2] = 3;
    fpd_sb_thrshd0_w0[3] = 3;
    fpd_sb_thrshd0_w0[4] = 3;
    fpd_sb_thrshd0_w0[5] = 3;
    fpd_sb_thrshd0_w0[6] = 3;
    fpd_sb_thrshd0_w0[7] = 3;
    for (int i=8; i<fdec_max_itr;i++)
        fpd_sb_thrshd0_w0[i] = 2;

    fpd_sb_thrshd0_w1[0] = 6;
    fpd_sb_thrshd0_w1[1] = 6;
    fpd_sb_thrshd0_w1[2] = 6;
    fpd_sb_thrshd0_w1[3] = 6;
    fpd_sb_thrshd0_w1[4] = 6;
    fpd_sb_thrshd0_w1[5] = 6;
    fpd_sb_thrshd0_w1[6] = 6;
    fpd_sb_thrshd0_w1[7] = 6;
    for (int i=8; i<fdec_max_itr;i++)
        fpd_sb_thrshd0_w1[i] = 6;

    fpd_sb_thrshd1_s0[0] = 0;
    fpd_sb_thrshd1_s0[1] = 0;
    fpd_sb_thrshd1_s0[2] = 0;
    fpd_sb_thrshd1_s0[3] = 0;
    fpd_sb_thrshd1_s0[4] = 0;
    fpd_sb_thrshd1_s0[5] = 0;
    fpd_sb_thrshd1_s0[6] = 0;
    fpd_sb_thrshd1_s0[7] = 0;
    for (int i=8; i<fdec_max_itr;i++)
        fpd_sb_thrshd1_s0[i] = 0;

    fpd_sb_thrshd1_s1[0] = 6;
    fpd_sb_thrshd1_s1[1] = 6;
    fpd_sb_thrshd1_s1[2] = 6;
    fpd_sb_thrshd1_s1[3] = 6;
    fpd_sb_thrshd1_s1[4] = 6;
    fpd_sb_thrshd1_s1[5] = 6;
    fpd_sb_thrshd1_s1[6] = 6;
    fpd_sb_thrshd1_s1[7] = 6;
    for (int i=8; i<fdec_max_itr;i++)
        fpd_sb_thrshd1_s1[i] = 6;

    fpd_sb_thrshd1_w0[0] = 0;
    fpd_sb_thrshd1_w0[1] = 0;
    fpd_sb_thrshd1_w0[2] = 0;
    fpd_sb_thrshd1_w0[3] = 0;
    fpd_sb_thrshd1_w0[4] = 0;
    fpd_sb_thrshd1_w0[5] = 0;
    fpd_sb_thrshd1_w0[6] = 0;
    fpd_sb_thrshd1_w0[7] = 0;
    for (int i=8; i<fdec_max_itr;i++)
        fpd_sb_thrshd1_w0[i] = 0;

    fpd_sb_thrshd1_w1[0] = 5;
    fpd_sb_thrshd1_w1[1] = 5;
    fpd_sb_thrshd1_w1[2] = 5;
    fpd_sb_thrshd1_w1[3] = 4;
    fpd_sb_thrshd1_w1[4] = 5;
    fpd_sb_thrshd1_w1[5] = 5;
    fpd_sb_thrshd1_w1[6] = 5;
    fpd_sb_thrshd1_w1[7] = 5;
    fpd_sb_thrshd1_w1[8] = 5;
    fpd_sb_thrshd1_w1[9] = 5;
   for (int i=8; i<fdec_max_itr;i++)
        fpd_sb_thrshd1_w1[i] = 5;

} // ldpc_dec_config

void ldpc_packet::ldpc_clean()
{
    mod2sparse_free(qc_bm);
    mod2sparse_free(qc_hm);
    mod2sparse_free(qc_a);
    mod2sparse_free(qc_b);
    mod2sparse_free(qc_c);
    mod2sparse_free(qc_d);
    mod2sparse_free(qc_e);
    mod2sparse_free(qc_fi);
    mod2sparse_free(qc_g);
    mod2sparse_free(qc_f1);
    mod2sparse_free(qc_f2);

    free(flp_thrshd0);
    flp_thrshd0 = NULL;
    free(flp_thrshd1);
    flp_thrshd1 = NULL;
    free(flp_thrshd0_s);
    flp_thrshd0_s = NULL;
    free(flp_thrshd1_s);
    flp_thrshd1_s = NULL;
    free(flp_thrshd0_w);
    flp_thrshd0_w = NULL;
    free(flp_thrshd1_w);
    flp_thrshd1_w = NULL;
    free(sb_thrshd0_s0);
    sb_thrshd0_s0 = NULL;
    free(sb_thrshd0_s1);
    sb_thrshd0_s1 = NULL;
    free(sb_thrshd1_s0);
    sb_thrshd1_s0 = NULL;
    free(sb_thrshd1_s1);
    sb_thrshd1_s1 = NULL;
    free(sb_thrshd0_w0);
    sb_thrshd0_w0 = NULL;
    free(sb_thrshd0_w1);
    sb_thrshd0_w1 = NULL;
    free(sb_thrshd1_w0);
    sb_thrshd1_w0 = NULL;
    free(sb_thrshd1_w1);
    sb_thrshd1_w1 = NULL;

    for (int i=0; i<bm_m; i++)
    {
        delete [] drop_col_bit_map[i];
        delete [] drop_col[i];
    }
    delete [] drop_col_bit_map;
    delete [] drop_col;
    delete [] wit_drop;
}

void ldpc_packet::ldpc_pckt_alloc()
{
    ch_pckt_alloc();

    dec_blk = (char *)calloc(blk_len, sizeof(*dec_blk));
    dec_do_blk = (char *)calloc(hm_n, sizeof(*dec_do_blk));
    usr_blk = (char *)calloc(info_len, sizeof(*usr_blk));
    enc_di_blk = (char *)calloc(hm_k, sizeof(*enc_di_blk));
    enc_do_blk = (char *)calloc(hm_n, sizeof(*enc_do_blk));
    dec_di_blk = (char *)calloc(hm_n, sizeof(*dec_di_blk));

    cw_fail = 0;
    cw_miscorr = 0;
    cor_err_num = 0;
    dec_err_num = 0;
}

void ldpc_packet::ldpc_pckt_clean()
{
    ch_pckt_clean();

    free(usr_blk);
    usr_blk = NULL;
    free(dec_blk);
    dec_blk = NULL;
    free(enc_di_blk);
    enc_di_blk = NULL;
    free(enc_do_blk);
    enc_do_blk = NULL;
    free(dec_di_blk);
    dec_di_blk = NULL;
    free(dec_do_blk);
    dec_do_blk = NULL;

} // ldpc_pckt_free

void ldpc_packet::ldpc_encoder()
{

}

void ldpc_packet::ldpc_decoder(enum dec_model dec_mode)
{
    // add 0 padding
    vec_copy(det_blk, dec_di_blk, 0, 0, info_len);
    for (int i=0; i<pad_len; i++)
        dec_di_blk[info_len + i] = max_llr_bin;

    vec_copy(det_blk, dec_di_blk, info_len, hm_k, hm_m);
    for (int i=hm_n-1;i>hm_n-cir_sz+pad_bit; i--)
        dec_di_blk[i] = 0;

    if (dec_mode == SKIP)
        ldpc_dec_skip();
    else if (dec_mode == BF_P0)
        ldpc_dec_bf(0, col_skip_itr);
    else if (dec_mode == BF_P3)
        ldpc_dec_bf(3, col_skip_itr);
    else if (dec_mode == BF_G2)
    {
        if (reg_fp_flg)
            ldpc_dec_bf3(3, col_skip_itr);
        else
            ldpc_dec_bf2(3, col_skip_itr);
    }
    else if (dec_mode == LAYER)
        ldpc_dec_layer();
    else if (dec_mode == TBFDEC)
        ldpc_dec_2bit_bf(3, col_skip_itr);

    //remove padding
    vec_copy(dec_do_blk, dec_blk, 0, 0, info_len);
    vec_copy(dec_do_blk, dec_blk, hm_k, info_len, hm_m);

    // check error bit number
    dec_err_num = 0;
    cor_err_num = 0;
    for (int i=0; i<real_len; i++)
    {
        if (dec_blk[i] != tx_blk[i])
            dec_err_num++;
        if (dec_blk[i] != det_blk[i])
            cor_err_num++;
    }

    // mis-correction case
    if ((cw_fail==0) && (dec_err_num!=0))
        cw_miscorr = 1;
    else
        cw_miscorr = 0;

#ifdef _LDPC_DEBUG
    if (cw_fail == 0)
    {
        if (dec_err_num == 0)
            printf("[LDPC DEBUG] Decoding success, no error bits found.\n");
        else
            printf("[LDPC DEBUG] MIS-CORRECTION!!! (%d)\n", dec_err_num);
    }
    else
    {
        printf("[LDPC DEBUG] Decoding failed, (%d)!\n", dec_err_num);
    }
#endif
}


void ldpc_packet::ldpc_dec_bf(int p_num, int col_skip_itr)
{
    mod2entry *e;
    int synd_wt;
    int col_updt;
    int itr_updt;
    bool col_skip;
    char *cn_synd_mem;  //syndrome memory in CN order
    char *cn_synd_sel;  //selected syndrome in CN order
    char *vn_synd_sel;  //selected syndrome in VN order
    char *vn_synd_cnt;  //syndrome weight of select columns in VN order
    char *vn_hd_sel;    //current HD of selected column (from dec_do_blk) in VN order
    char *vn_raw_sel;    //raw data of selected column (from dec_di_blk) in VN order
    char **vn_flp_sel;  //flip flag of selected column in VN order
    int *vn_flp_col;    //column index of the vn_flp_sel
    int *vn_flp_itr;    //iteration of the vn_flp_sel
    char *cn_flp_sel;  //flip flag of selected column in CN order
    char *cn_synd_new;   // new syndrome in CN order
    char *cn_synd_old;  


    //allocate memory
    cn_synd_mem = (char *)calloc(hm_m, sizeof(*cn_synd_mem));
    cn_synd_sel = (char *)calloc(cir_sz, sizeof(*cn_synd_sel));
    vn_synd_sel = (char *)calloc(cir_sz, sizeof(*vn_synd_sel));
    vn_synd_cnt = (char *)calloc(cir_sz, sizeof(*vn_synd_cnt));
    vn_hd_sel = (char *)calloc(cir_sz, sizeof(*vn_hd_sel));
    vn_raw_sel = (char *)calloc(cir_sz, sizeof(*vn_raw_sel));
    vn_flp_sel = (char **)calloc(p_num+1, sizeof(*vn_flp_sel));
    for(int i=0;i<=p_num;i++)
        vn_flp_sel[i] = (char *)calloc(cir_sz, sizeof(*vn_flp_sel[i]));
    vn_flp_col = (int *)calloc(p_num+1, sizeof(*vn_flp_col));
    vn_flp_itr = (int *)calloc(p_num+1, sizeof(*vn_flp_itr));
    cn_flp_sel = (char *)calloc(cir_sz, sizeof(*cn_flp_sel));
    cn_synd_new = (char *)calloc(cir_sz, sizeof(*cn_synd_new));
    cn_synd_old = (char *)calloc(cir_sz, sizeof(*cn_synd_old));


    cw_fail = 1;
    cw_miscorr = 0;
    fina_synd_wt =0;


    vec_copy(dec_di_blk, dec_do_blk, 0,0,hm_n);
    vec_clr(cn_synd_mem,hm_m);
    for(int i=0;i<p_num;i++)
        vec_clr(vn_flp_sel[i], cir_sz);
    fdec_cyc_num = 0;
    fdec_cyc_org = 0;


    for(int itr=0;(itr<=fdec_max_itr)&&((fdec_early_term_en==0)||(cw_fail==1));itr++)
    {
        for(int i=0;i<(itr==fdec_max_itr?p_num:bm_n);i++)
        {
            if(itr==0)
            {
                vec_copy(dec_di_blk,vn_flp_sel[p_num],i*cir_sz,0,cir_sz);
                col_updt=i;
                itr_updt = itr;
                col_skip=false;
            }
            else
            {
                vec_clr(vn_synd_cnt,cir_sz);
                for(e=mod2sparse_first_in_col(qc_bm,i);!mod2sparse_at_end(e);e=mod2sparse_next_in_col(e))
                {
                    int drop_start = cir_sz +1;
                    int drop_length = 0;
                    int min_drop_ind = -1;
                    vec_copy(cn_synd_mem,cn_synd_sel,e->row*cir_sz,0,cir_sz);


                    if(e->col==bm_n-1)
                    {
                        drop_start = (cir_sz - (e->shift) + pad_bit) % cir_sz;
                        drop_length = cir_sz - pad_bit;
                        min_drop_ind = drop_start +drop_length;
                    }
                    else
                    {
                        if(drop_col_bit_map[e->row][e->col])
                        {
                            if(e->row == bm_m -1)
                            {
                                drop_start = pad_bit;
                                drop_length = cir_sz - pad_bit;
                            }
                            else
                            {
                                drop_start = (e->col - e->row)/(bm_m-1)*pad_bit % cir_sz;
                                drop_length = pad_bit;
                            }
                            min_drop_ind = drop_start+drop_length;
                        }
                    }


                    for(int ii=0;ii<cir_sz;ii++)
                    {
                        if((min_drop_ind<=cir_sz && ii>=drop_start && ii < min_drop_ind)
                        || (min_drop_ind>cir_sz && (ii>=drop_start || ii<min_drop_ind % cir_sz)))
                            cn_synd_sel[ii] = 0;
                    }


                    vec_shift(cn_synd_sel,vn_hd_sel,cir_sz,e->shift);
                    vec_incr(vn_synd_cnt, vn_synd_sel,cir_sz);
                }


                // previous column is skipped
                if(col_skip)
                    col_skip = false;   // Column skip feature OFF
                else
                    if(col_skip_itr==0)
                        col_skip = false;// non-skip iterations
                    else
                        if((col_skip_itr>0)&&(itr<col_skip))
                            col_skip = false;
                        else
                        {
                            col_skip = true;
                            //make a skip decision
                            for(int j=0;j<cir_sz;j++)
                            {
                                if(vn_synd_cnt[j]>=flp_thrshd1[itr-1])
                                    col_skip = false;
                            }
                        }
               
                //flip logics
                if(col_skip == false)
                {
                    // read raw and current HD
                    vec_copy(dec_di_blk,vn_raw_sel, i*cir_sz,0,cir_sz);
                    vec_copy(dec_do_blk,vn_hd_sel, i*cir_sz,0,cir_sz);


                    //pipelines
                    for(int j=p_num;j>0;j--)
                        vec_copy(vn_flp_sel[j-1],vn_flp_sel[j],0,0,cir_sz);
                    for(int j=p_num;j>0;j--)
                    {
                        vn_flp_col[j] = vn_flp_col[j-1];
                        vn_flp_itr[j] = vn_flp_itr[j-1];
                    }
                    vn_flp_col[0] = i;
                    vn_flp_itr[0] = itr;


                    for(int j=0;j<cir_sz;j++)
                    {
                        if(((vn_raw_sel[j]==vn_hd_sel[j]) && (vn_synd_cnt[j]>=flp_thrshd0[itr-1]))
                        || ((vn_raw_sel[j]!=vn_hd_sel[j]) && (vn_synd_cnt[j]>=flp_thrshd1[itr-1])))
                        {
                            vn_flp_sel[0][j] = 1;
                        }
                        else
                        {
                            vn_flp_sel[0][j] = 0;
                        }
                    }


                    col_updt = vn_flp_col[p_num];
                    itr_updt = vn_flp_itr[p_num];


                    for(int j=0;j<cir_sz;j++)
                    {
                        if(vn_flp_sel[p_num][j] == 1)
                        {
                            dec_do_blk[col_updt*cir_sz+j] = (dec_do_blk[col_updt*cir_sz+j]+1)%2;
                        }
                    }
                }//non-skipped columns(flip logic)
            }// non-1st iteration columns
            fdec_cyc_org++;


            if(col_skip==false)
            {
                fdec_cyc_num++;


                for(e=mod2sparse_first_in_col(qc_bm,col_updt);!mod2sparse_at_end(e);e=mod2sparse_next_in_col(e))
                {
                    int drop_start = cir_sz+1;
                    int drop_length = 0;
                    int min_drop_ind =-1;


                    vec_shift(vn_flp_sel[p_num], cn_flp_sel,cir_sz,-1*e->shift);
                    if(col_updt==bm_n-1)
                    {
                        drop_start = (cir_sz - (e->shift) + pad_bit) % cir_sz;
                        drop_length = cir_sz - pad_bit;
                        min_drop_ind = drop_start + drop_length;
                    }
                    else
                    {
                        if(drop_col_bit_map[e->row][e->col])
                        {
                            if(e->row==bm_m-1)
                            {
                                drop_start = pad_bit;
                                drop_length = cir_sz - pad_bit;
                            }
                            else
                            {
                                drop_start = (e->col-e->row)/(bm_m-1) * pad_bit % cir_sz;
                                drop_length = pad_bit;
                            }
                            min_drop_ind = drop_start + drop_length;
                        }
                    }
                    for(int ii=0;ii<cir_sz;ii++)
                    {
                        if((min_drop_ind<=cir_sz && ii>=drop_start && ii < min_drop_ind)
                        || (min_drop_ind>cir_sz && (ii>=drop_start || ii<min_drop_ind % cir_sz)))
                            cn_flp_sel[ii] = 0;
                    }
                    vec_copy(cn_synd_mem,cn_synd_old,e->row*cir_sz,0,cir_sz);
                    vec_mod2_add(cn_flp_sel, cn_synd_old, cn_synd_new, cir_sz);
                    vec_copy(cn_synd_new,cn_synd_mem,0,e->row*cir_sz,cir_sz);


                }


                if((itr>0) || (i==(bm_n-1)))
                {
                    synd_wt = vec_sum(cn_synd_mem,hm_m);
                    if(synd_wt ==0)
                    {
                        cw_fail=1;
                        cnvg_itr = itr_updt;
                        cnvg_lyr = col_updt;
                        if(fdec_early_term_en ==1)
                            break;
                    }
                }
            }
        }
    }


    if((fdec_early_term_en == 0) || (cw_fail==1))
    {
        cnvg_itr = fdec_max_itr - 1;
        cnvg_lyr = bm_n - 1;
        fina_synd_wt = vec_sum(cn_synd_mem,hm_m);
    }


    free(cn_synd_mem);
    free(cn_synd_sel);
    free(vn_synd_sel);
    free(vn_synd_cnt);
    free(vn_hd_sel);
    free(vn_raw_sel);
    for(int i=0;i<=p_num;i++)
        free(vn_flp_sel[i]);
    free(vn_flp_sel);
    free(vn_flp_col);
    free(vn_flp_itr);
    free(cn_flp_sel);
    free(cn_synd_new);
    free(cn_synd_old);

#ifdef _LDPC_DEBUG_DUMP
    fclose(cfp);
    fclose(sfp);
    fclose(hdfp);
    fclose(lfp);
#endif
}


void ldpc_packet::ldpc_dec_layer()
{
#ifdef _LDPC_DEBUG_DUMP
    FILE *cfp, *sfp, *hdfp, *lfp;
    int stmp, vtmp;
    char cmem_dump[50] = "./output/rdec_cmem_dump.txt";
    char stot_dump[50] = "./output/rdec_stot_dump.txt";
    char hdmem_dump[50] = "./output/rdec_hdmem_dump.txt";
    char log_dump[50] = "./output/rdec_log_dump.txt";
    cfp = fopen(cmem_dump, "w");
    sfp = fopen(stot_dump, "w");
    hdfp = fopen(hdmem_dump, "w");
    lfp = fopen(log_dump, "w");
#endif

    mod2entry *e, *e_pre;
    char *dec_init;
    int shift_val1;
    int shift_val2;
    int cir_cnt;
    int sign_tmp;
    float val_tmp;
    int hd_init;

    struct cn_msg **cn_c_mem;
    struct cn_msg *cn_c_updt_cur;   // current layer check node msg to be updt
    struct cn_msg *cn_c_sel_cur;    // current layer check node msg
    struct cn_msg *cn_c_sel_pre;   // previous layer check node msg
    float **cn_q_mem;       // Q mem in CN order of previous layer
    float *cn_q_sel_pre;    // Q msg of the select circulant from previous layer
    float *cn_q_sel_cur;    // Q msg of the select circulant from current layer
    float *cn_r_new_pre;    // New R msg in CN order of previous layer
    float *cn_r_old_cur;    // old R msg in CN order of previous layer
    float *cn_app_pre;      // APP = Q + R_new in CN order of previous layer
    float *cn_app_cur;      // APP = Q + R_new in CN order of current layer
    float *cn_q_updt_cur;   // Updated Q msg of the select circulant in current layer
    int **cn_q_sign;      // Q sign

    char *layer_synd;
    char *cn_dec_hd;
    char *vn_dec_hd;
    int hd_updated;
    int layer_synd_wt;
    int synd_pass_cnt = 0;
    int hd_stable_cnt = 0;
    char *drop_wit_now;
    
    // allocation
    dec_init = (char *)calloc(hm_n, sizeof(*dec_init));
    drop_wit_now = (char *)calloc(bm_m, sizeof(*drop_wit_now));
    vec_set(dec_init, bm_n);

    cn_c_mem = (struct cn_msg **)calloc(bm_m, sizeof(*cn_c_mem));
    for (int i = 0; i<bm_m*col_wt; i++)
        cn_q_sign[i] = (int*)calloc(cir_sz, sizeof(*cn_q_sign[i]));

    layer_synd = (char *)calloc(cir_sz, sizeof(*layer_synd));
    vn_dec_hd = (char *)calloc(cir_sz, sizeof(*vn_dec_hd));
    cn_dec_hd = (char *)calloc(cir_sz, sizeof(*cn_dec_hd));

    // initialize decoder
    cw_fail = 1;
    cw_miscorr = 0;
    vec_copy(dec_di_blk, dec_do_blk, 0, 0, hm_n);

    for (int i = 0; i < bm_n; i++)
    {
        for (int j=0; j<cir_sz; j++)
        {
            cn_q_mem[i][j] = (float)llr_tbl[dec_di_blk[i*cir_sz+j]];

            // Quantization
            if (finite_mode == 1)
            {
                cn_q_mem[i][j] = (float)Sat_Quan((double)cn_q_mem[i][j], finite_q_max, finite_q_min, finite_q_num, finite_f_num);
            }
        }
    }

    // iterative decoding
    for(int itr=0;(itr<=ldec_max_itr)&&((ldec_early_term_en==0)||(cw_fail==1));itr++)
    {
        // Q sign mem index
        cir_cnt = 0;
        vec_clr(drop_wit_now, bm_m);

        // layer decoding
        for(int layer=0;layer<bm_m &&((ldec_early_term_en==0)||(cw_fail==1));layer++)
        {
            // initilize HD mem
            hd_init = (vec_sum(dec_init, bm_n)!=0);

#ifdef _LDPC_DEBUG_DUMP
            printf("[LDPC DEBUG] Layer decoding @ iteration %d, layer %d ...\n", itr, layer);
#endif                        

            // init current layer C-MSG
            // C-MSG of previous iteration

            cn_c_sel_cur = cn_c_mem[layer];
            // C-MSG to be updt
            for (int i = 0; i < cir_sz; i++)
            {
                cn_c_updt_cur[i].min1_val = 100000;
                cn_c_updt_cur[i].min2_val = 100000;
                cn_c_updt_cur[i].min1_pos = 0;
                cn_c_updt_cur[i].sign_tot = 1;
            }

            // Earlier termination init
            hd_updated = 0;
            vec_clr(layer_synd, cir_sz);

            // Per circulant of the layer
            for (e = mod2sparse_first_in_row(qc_bm, layer);
                !mod2sparse_at_end(e) &&((ldec_early_term_en==0)||(cw_fail==1));
                e = mod2sparse_next_in_row(e))
            {
                // read Q from the previous layer of the selected column
                cn_q_sel_pre = cn_q_mem[e->col];
                int drop_start = cir_sz + 1;
                int drop_length = 0;
                int mem_app_start = cir_sz + 1;
                int mem_app_length = 0;
                int min_app_ind = -1;
                int min_drop_ind = -1;

                if(drop_col[e->row][drop_wit_now[e->row]] == e->col && drop_wit_now[e->row] < wit_drop[e->row])
                {
                    if (e->row == bm_m - 1)
                    {
                        drop_start = pad_bit;
                        drop_length = cir_sz - pad_bit;
                    }
                    else
                    {
                        drop_start = drop_wit_now[e->row] * pad_bit % cir_sz;
                        drop_length = pad_bit;
                    }
                    min_drop_ind = drop_start + drop_length;
                    drop_wit_now[e->row]++;
                }

                if(e->col == bm_n - 1)
                {
                    drop_start = (cir_sz - (e->shift) + pad_bit) % cir_sz;
                    drop_length = cir_sz - pad_bit;
                    min_drop_ind = drop_start + drop_length;
                }

                // find out the previous layer of select column
                e_pre = mod2sparse_prev_in_col(e);
                if (mod2sparse_at_end(e_pre))
                    e_pre = mod2sparse_last_in_col(qc_bm, e->col);
                
                if ((e_pre->row == bm_m - 1 || e_pre->col%(bm_m-1)==e_pre->row)
                && (drop_col[e_pre->row][wit_drop[e_pre->row]-1] >= e_pre->col))
                {
                    if (e_pre->row == bm_m - 1)
                    {
                        mem_app_start = pad_bit;
                        mem_app_length = cir_sz - pad_bit;
                    }
                    else
                    {
                        mem_app_start = (e_pre->col - e_pre->row) / (bm_m - 1) * pad_bit % cir_sz;
                        mem_app_length = pad_bit;
                    }
                    min_app_ind = (mem_app_start + mem_app_length);
                }

                if (e_pre->col == bm_n - 1)
                {
                    mem_app_start = (cir_sz - (e_pre->shift) + pad_bit) % cir_sz;
                    mem_app_length = (cir_sz - pad_bit) % cir_sz;
                    min_app_ind = mem_app_start + mem_app_length;
                }

                cn_c_sel_pre = cn_c_mem[e_pre->row];    // read previous layer C msg

                // cal Rnew and APP
                for (int i = 0; i < cir_sz; i++)
                {
                    // Qmsg sign
                    sign_tmp = (cn_q_sel_pre[i] < 0) ? -1 : 1;
                    if ((min_app_ind <= cir_sz && i >= mem_app_start && i < min_app_ind)
                    || (min_app_ind > cir_sz && (i >= mem_app_start || i <min_app_ind % cir_sz)))
                        cn_app_pre[i] = 0;
                    else
                    {
                        if (cn_c_sel_pre[i].min1_pos == e->col)
                        {
                            cn_r_new_pre[i] = cn_c_sel_pre[i].min2_val * cn_c_sel_pre[i].sign_tot * sign_tmp;
                        }
                        else
                        {
                            cn_r_new_pre[i] = cn_c_sel_pre[i].min1_val * cn_c_sel_pre[i].sign_tot * sign_tmp;
                        }
                    }

                    // Rnew
                    // APP in CN order of previous layer
                    cn_app_pre[i] = cn_q_sel_pre[i] + cn_r_new_pre[i];

                    if (finite_mode == 1)
                    {
                        cn_app_pre[i] = (float)Sat_Quan((double)cn_app_pre[i], finite_q_max, finite_q_min, finite_q_num, finite_f_num);
                    }
                }

                // APP shift
                // when decoder initilized, Q msg are in VN order
                if (dec_init[e->col] == 1)
                {
                    shift_val1 = e->shift;
                    shift_val2 = 0;
                    dec_init[e->col] = 0;
                }
                else 
                {
                    shift_val1 = -1 * e_pre->shift + e->shift;
                    shift_val2 = -1 * e_pre->shift;
                }

                for (int i = 0; i < cir_sz; i++)
                {
                    cn_app_cur[i] = cn_app_pre[(i+shift_val1 + cir_sz) % cir_sz];
                    vn_dec_hd[i] = cn_app_pre[(i+shift_val2 + cir_sz) % cir_sz] > 0 ? 0 : 1;
                }

#ifdef _LDPC_DEBUG_DUMP
                fprintf(lfp, "ITR%2d/LAYER%2d/COL%2d: \n", itr, layer, e->col);
                for (int i=0; i<cir_sz; i++)
                {
                    fprintf(lfp, "Q PRE MSG:");

                    for (int j=0; j<8; j++)
                    {
                        vtmp = int(cn_q_sel_pre[i*8+j] * pow(2, finite_f_num));

                    if (vtmp < 0)
                    {
                        vtmp = -vtmp;
                        fprintf(lfp, " %3d-%2X", i*8+j, vtmp);
                    }
                    else
                        fprintf(lfp, " %3d+%2X", i*8+j, vtmp);
                    }
                    fprintf(lfp, "\n");
                }

                for (int i=0; i<cir_sz/8; i++)
                {
                    fprintf(lfp, "R NEW MSG:");

                    for (int j=0; j<8; j++)
                    {
                        vtmp = int(cn_r_new_pre[i*8+j] * pow(2, finite_f_num));

                        if (vtmp < 0)
                        {
                            vtmp = -vtmp;
                            fprintf(lfp, " %3d-%2X", i*8+j, vtmp);
                        }
                        else
                            fprintf(lfp, " %3d+%2X", i*8+j, vtmp);
                    }
                    fprintf(lfp, "\n");
                }

                for (int i=0; i<cir_sz/8; i++)
                {
                    fprintf(lfp, "APP-C MSG:");

                    for (int j=0; j<8; j++)
                    {
                        vtmp = int(cn_app_pre[i*8+j] * pow(2, finite_f_num));

                        if (vtmp < 0)
                        {
                            vtmp = -vtmp;
                            fprintf(lfp, " %3d-%2X", i*8+j, vtmp);
                        }
                        else
                            fprintf(lfp, " %3d+%2X", i*8+j, vtmp);
                    }
                    fprintf(lfp, "\n");
                }

                for (int i=0; i<cir_sz/8; i++)
                {
                    fprintf(lfp, "APP-S MSG:");

                    for (int j=0; j<8; j++)
                    {
                        vtmp = int(cn_app_cur[i*8+j] * pow(2, finite_f_num));

                        if (vtmp < 0)
                        {
                            vtmp = -vtmp;
                            fprintf(lfp, " %3d-%2X", i*8+j, vtmp);
                        }
                        else
                            fprintf(lfp, " %3d+%2X", i*8+j, vtmp);
                    }
                    fprintf(lfp, "\n");
                }
#endif
                // CW converge check logic per circulant 
                // 1. check if HD updated
                if (hd_updated == 0)
                {
                    if (vec_cmp(dec_do_blk, vn_dec_hd, e->col*cir_sz, 0, cir_sz) == 1)
                    {
                        hd_updated = 1;
                    }
                }
#ifdef _LDPC_DEBUG_DUMP
                for (int i = 0; i < cir_sz; i++)
                {
                    if (dec_do_blk[e->col*cir_sz + i] != vn_dec_hd[i])
                    {
                        fprintf(lfp, "ITR%2d/LAYER%2d/COL%wd: flip bit %d (%d-->%d)\n", itr, layer, e->col, i, dec_do_blk[e->col*cir_sz + i], vn_dec_hd[i]);
                    }
                }
#endif
                vec_copy(vn_dec_hd, dec_do_blk, 0, e->col*cir_sz, cir_sz);

#ifdef _LDPC_DEBUG_DUMP
                fprintf(hdfp, "ITR%2d/LAYER%2d/COL%wd: ", itr, layer, e->col);
                for (int i=cir_sz/4-1; i>=0; i--)
                {
                    stmp = 0;
                    for (int j=3; j>=0; j--)
                    {
                        stmp = stmp *2 +vn_dec_hd[i*4+j];
                    }
                    fprintf(hdfp, "%lx", stmp);
                }
                fprintf(hdfp, "\n");
#endif

                // 2. accumulate syndrome
                vec_shift(vn_dec_hd, cn_dec_hd, cir_sz, -1 * e->shift);

                for (int i = 0; i < cir_sz; i++)
                {
                    if ((min_drop_ind <= cir_sz && i >= drop_start && i < min_drop_ind)
                    || (min_drop_ind > cir_sz && (i >= drop_start || i < min_drop_ind % cir_sz)))
                        cn_dec_hd[i] = 0;
                }
                vec_mod2_add(cn_dec_hd, layer_synd, layer_synd, cir_sz);

                // calculate R_old and current Q, update current layer C and Q
                for (int i=0; i<cir_sz; i++)
                {
                    if ((min_drop_ind <= cir_sz && i >= drop_start && i < min_drop_ind)
                    || (min_drop_ind > cir_sz && (i >= drop_start || i < min_drop_ind % cir_sz)))
                        cn_r_old_cur[i] = 0;
                    else
                    {
                        if (cn_c_sel_cur[i].min1_pos == e->col)
                        {
                            cn_r_old_cur[i] = cn_c_sel_cur[i].min2_val * cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];
                        }
                        else
                        {
                            cn_r_old_cur[i] = cn_c_sel_cur[i].min1_val * cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];
                        }
                    }

                    // Q -= Rold
                    cn_q_updt_cur[i] = cn_app_cur[i] - cn_r_old_cur[i];

                    if (finite_mode == 1)
                    {
                        cn_q_updt_cur[i] = (float)Sat_Quan((double)cn_q_updt_cur[i], finite_q_max, finite_q_min, finite_q_num, finite_f_num);
                    }

                    // update C

                    if ((min_drop_ind <= cir_sz && i >= drop_start && i < min_drop_ind)
                    || (min_drop_ind > cir_sz && (i >= drop_start || i < min_drop_ind % cir_sz)))
                    {
                    }
                    else
                    {
                        sign_tmp = (cn_q_updt_cur[i] < 0) ? -1 : 1;
                        val_tmp = cn_q_updt_cur[i] * sign_tmp;
                        cn_c_updt_cur[i].sign_tot *= sign_tmp;

                        if (val_tmp < cn_c_updt_cur[i].min1_val)
                        {
                            cn_c_updt_cur[i].min2_val = cn_c_updt_cur[i].min1_val;
                            cn_c_updt_cur[i].min1_val = val_tmp;
                            cn_c_updt_cur[i].min1_pos = e->col;
                        }
                        else if (val_tmp < cn_c_updt_cur[i].min2_val)
                        {
                            cn_c_updt_cur[i].min2_val = val_tmp;
                        }

                        // update Q sign
                        cn_q_sign[cir_cnt][i] = sign_tmp;
                    }
                }

                // update Q memory
                for (int i=0; i<cir_sz; i++)
                {
                    cn_q_mem[e->col][i] = cn_q_updt_cur[i];
                }

#ifdef _LDPC_DEBUG_DUMP
                for (int i=0; i<cir_sz/8; i++)
                {
                    fprintf(lfp, "R OLD MSG:");

                    for (int j=0; j<8; j++)
                    {
                        vtmp = int(cn_r_old_cur[i*8+j] * pow(2, finite_f_num));

                        if (vtmp < 0)
                        {
                            vtmp = -vtmp;
                            fprintf(lfp, " %3d-%2X", i*8+j, vtmp);
                        }
                        else
                            fprintf(lfp, " %3d+%2X", i*8+j, vtmp);
                    }
                    fprintf(lfp, "\n");
                }

                for (int i=0; i<cir_sz/8; i++)
                {
                    fprintf(lfp, "Q NEW MSG:");

                    for (int j=0; j<8; j++)
                    {
                        vtmp = int(cn_q_updt_cur[i*8+j] * pow(2, finite_f_num));

                        if (vtmp < 0)
                        {
                            vtmp = -vtmp;
                            fprintf(lfp, " %3d-%2X", i*8+j, vtmp);
                        }
                        else
                            fprintf(lfp, " %3d+%2X", i*8+j, vtmp);
                    }
                    fprintf(lfp, "\n");
                }
#endif
                cir_cnt++;
            } // per circulant

            // update C_MSG per layer
            for (int i=0; i<cir_sz; i++)
            {
                cn_c_mem[layer][i].min1_val = cn_c_updt_cur[i].min1_val * alpha;
                cn_c_mem[layer][i].min2_val = cn_c_updt_cur[i].min2_val * alpha;

                if (finite_mode == 1)
                {
                    cn_c_mem[layer][i].min1_val = (float)Sat_Quan((double)cn_c_mem[layer][i].min1_val, finite_q_max, finite_q_min, finite_q_num, finite_f_num);
                    cn_c_mem[layer][i].min2_val = (float)Sat_Quan((double)cn_c_mem[layer][i].min2_val, finite_q_max, finite_q_min, finite_q_num, finite_f_num);
                }

                cn_c_mem[layer][i].min1_pos = cn_c_updt_cur[i].min1_pos;
                cn_c_mem[layer][i].sign_tot = cn_c_updt_cur[i].sign_tot;

#ifdef _LDPC_DEBUG_DUMP
                fprintf(cfp, "ITR%d/L%d/C%d: min1 %x, min2 %x, min1 pos %d, sign_tot %d\n",
                itr, layer, i, int(16*cn_c_mem[layer][i].min1_val), int(16*cn_c_mem[layer][i].min2_val), cn_c_mem[layer][i].min1_pos, (1-cn_c_mem[layer][i].sign_tot)/2);
#endif
            }

#ifdef _LDPC_DEBUG_DUMP
            for (int i=(cir_sz/4-1); i>=0; i--)
            {
                stmp = 0;
                for (int j=3; j>=0; j--)
                {
                    stmp = stmp *2 + (1-cn_c_mem[layer][i*4+j].sign_tot)/2;
                }
                fprintf(sfp, "%lx", stmp);
            }
            fprintf(sfp, "\n");
#endif

            // check converage checking
            layer_synd_wt = vec_sum(layer_synd, cir_sz);
            if (hd_init == 1)
            {
                hd_stable_cnt = 0;
                synd_pass_cnt = 0;

#ifdef _LDPC_DEBUG
                printf("[LDPC DEBUG] HD stable and syndrome passed @ iteration %d, layer %d\n", itr, layer);
#endif
#ifdef _LDPC_DEBUG_DUMP
                fprintf(lfp, "ITR%d/LAYER%d: HD stable for %d and syndrome passed for %d layers\n", itr, layer, hd_stable_cnt, synd_pass_cnt);
#endif
                hd_stable_cnt++;
                synd_pass_cnt++;
            }
            else
            {
                hd_stable_cnt = 0;
                synd_pass_cnt = 0;
#ifdef _LDPC_DEBUG_DUMP
                if (hd_updated == 1)
                    fprintf(lfp, "ITR%d/LAYER%d: HD memory is updated\n", itr, layer);

                if (layer_synd_wt == 1)
                    fprintf(lfp, "ITR%d/LAYER%d: Syndrome check fail\n", itr, layer);
#endif
            }
#ifdef _LDPC_DEBUG
                printf("[LDPC DEBUG] HD stable and syndrome passed @ iteration %d, layer %d\n", hd_stable_cnt, synd_pass_cnt);
#endif

            if ((hd_stable_cnt >= bm_m -1) && (synd_pass_cnt >= bm_m))
            {
                cw_fail = 0;
                cnvg_itr = itr;
                cnvg_lyr = layer;
#ifdef _LDPC_DEBUG
                printf("[LDPC DEBUG] Layer decoding converged @ iteration %d, layer %d\n", itr, layer);
#endif
            }
        }
    }

    if ((cw_fail == 0) || (ldec_early_term_en == 0))
    {
        cnvg_itr = ldec_max_itr - 1;
        cnvg_lyr = bm_m - 1;
    }

    // free all
    free(dec_init);
    for (int i=0; i<bm_m; i++)
        free(cn_c_mem[i]);
    free(cn_c_mem);
    free(cn_c_updt_cur);
    for (int i=0; i<bm_m*col_wt; i++)
        free(cn_q_sign[i]);
    for (int i=0; i<bm_n; i++)
        free(cn_q_mem[i]);
    free(cn_q_mem);
    free(cn_r_new_pre);
    free(cn_r_old_cur);
    free(cn_app_pre);
    free(cn_app_cur);
    free(cn_q_sel_cur);
    free(cn_q_sel_pre);
    free(cn_q_updt_cur);
    free(cn_dec_hd);
    free(vn_dec_hd);
    free(layer_synd);
    free(drop_wit_now);

#ifdef _LDPC_DEBUG_DUMP
    fclose(cfp);
    fclose(sfp);
    fclose(hdfp);
    fclose(lfp);
#endif
} // ldpc_dec_layer


void ldpc_packet::ldpc_dec_skip()
{
    for (int i=0; i<hm_n; i++)
        dec_do_blk[i] = dec_di_blk[i];

    cw_fail = ldpc_synd(dec_do_blk);
}

int ldpc_packet::ldpc_synd(char *cw)
{
    int synd_fail = 0;
    char *synd;

    synd = (char *)calloc(hm_m, sizeof(*synd));

    mod2sparse_mulvec(qc_hm, cw, synd);

    for (int i=0; i<hm_m; i++)
    {
        if (synd[i] != 0)
        {
            synd_fail = 1;
            break;
        }
    }

    free(synd);

    return synd_fail;
}