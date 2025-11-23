#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "finite_lib.h"
#include "ldpc_codec.h"
#include "mod2convert.h"
#include "mod2dense.h"
#include "mod2sparse.h"
#include "vec_op.h"

typedef struct
{
    int idx;
    int energy;
} mbf_candidate_t;

static int compare_mbf_candidate_desc(const void *a, const void *b)
{
    const mbf_candidate_t *ca = (const mbf_candidate_t *)a;
    const mbf_candidate_t *cb = (const mbf_candidate_t *)b;

    if (cb->energy != ca->energy)
        return (cb->energy - ca->energy);
    return (ca->idx - cb->idx);
}

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
        exit(0);
    }


    drop_col = new int *[bm_m];
    for(int i=0;i<bm_m-1;i++)
        drop_col[i] = new int [8];
    drop_col_bit_map = new int*[bm_m];
    drop_col[bm_m-1] = new int [80];
    for(int i=0;i<bm_m;i++)
    {
        drop_col_bit_map[i] = new int[bm_n];
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
                if(((i!=(j+tm_sz-bm_n+1))&&(col_shift>=0)) || ((i==(j+tm_sz-bm_n+1)) &&(col_shift!=0)))
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
                        mod2sparse_insert(qc_hm, i*cir_sz+k,j*cir_sz+(k+col_shift) % cir_sz);
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
        qc_t_cols[i] = i + (bm_n-tm_sz-1)*cir_sz;
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

#ifdef _LDPC_DUMP
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
    mod2dense_to_sparse(qc_fi_d, qc_fi);

    printf("[LDPC] Inverse F matrix generated!\n");

#ifdef _LDPC_FI_DUMP
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

#ifdef _LDPC_DUMP
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

#ifdef _LDPC_DEBUG
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

    tbbf_thrshd0 = (int *)calloc(128, sizeof(*flp_thrshd0));
    tbbf_thrshd1 = (int *)calloc(128, sizeof(*flp_thrshd0));
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
        flp_thrshd0_w[i] = 4;

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
        sb_thrshd0_w0[i] = 2;

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

} // ldpc_pckt_free

void ldpc_packet::ldpc_encoder()
{
    char *az1, *eaz1, *cz1, *sumz1, *z2, *bz2, *z3, *f2z, *f1z, *gz;

    az1 = (char *)calloc(tm_sz*cir_sz, sizeof(*az1));;
    eaz1 = (char *)calloc((bm_m-tm_sz-1)*cir_sz, sizeof(*eaz1));
    cz1 = (char *)calloc((bm_m-tm_sz-1)*cir_sz, sizeof(*cz1));
    sumz1 = (char *)calloc((bm_m-tm_sz-1)*cir_sz, sizeof(*sumz1));
    z2 = (char *)calloc((bm_m-tm_sz-1)*cir_sz, sizeof(*z2));
    bz2 = (char *)calloc(tm_sz*cir_sz, sizeof(*bz2));
    z3 = (char *)calloc(tm_sz*cir_sz, sizeof(*z3));
    f1z = (char *)calloc(tm_sz*cir_sz, sizeof(*f1z));
    f2z = (char *)calloc((bm_m-tm_sz-1)*cir_sz, sizeof(*f2z));
    gz = (char *)calloc(cir_sz, sizeof(*gz));

    // padding 0s
    vec_copy(usr_blk, enc_di_blk, 0, 0, info_len);
    for (int i=0; i<pad_len; i++)
        enc_di_blk[hm_k-pad_len+i] = 0;

    mod2entry *e;

    // A*Z1
    mod2sparse_mulvec(qc_a, enc_di_blk, az1);
    char *vn_flp_sel = (char *)calloc(cir_sz, sizeof(*vn_flp_sel));
    char *cn_flp_sel = (char *)calloc(cir_sz, sizeof(*cn_flp_sel));
    char *cn_synd_mem = (char *)calloc(hm_m, sizeof(*cn_synd_mem));
    char *cn_synd_old = (char *)calloc(cir_sz, sizeof(*cn_synd_old));
    char *cn_synd_new = (char *)calloc(cir_sz, sizeof(*cn_synd_new));
    vec_clr(cn_synd_mem, hm_m);
    for (int i=0; i<bm_k; i++)
    {
        vec_copy(enc_di_blk, vn_flp_sel, i*cir_sz, 0, cir_sz);
        for (e=mod2sparse_first_in_col(qc_bm, i); !mod2sparse_at_end(e); e=mod2sparse_next_in_col(e))
        {
            int drop_start = cir_sz + 1;
            int drop_length = 0;
            int min_drop_ind = -1;

            // barrel shifter
            vec_shift(vn_flp_sel, cn_flp_sel, cir_sz, -1 * e->shift);
            // read old syndrome
            if (i==bm_n-1)
            {
                drop_start = (cir_sz - (e->shift) + pad_bit) % cir_sz;
                drop_length = cir_sz - pad_bit;
                min_drop_ind = drop_start + drop_length;
            }
            else 
            {
                if (drop_col_bit_map[e->row][e->col] == 1)
                {
                    if (e->row == bm_m - 1)
                    {
                        drop_start = pad_bit;
                        drop_length = cir_sz - pad_bit;
                    }
                    else
                    {
                        drop_start = (e->col - e->row) / (bm_m - 1) * pad_bit %cir_sz;
                        drop_length = pad_bit;
                    }
                    min_drop_ind = drop_start + drop_length;
                }
            }
            for (int ii = 0; ii < cir_sz; ii++)
            {
                if ((min_drop_ind <= cir_sz && ii >= drop_start && ii < min_drop_ind) ||
                    (min_drop_ind >= cir_sz && (ii >= drop_start || ii < min_drop_ind % cir_sz)))
                    cn_flp_sel[ii] = 0;
            }

            vec_copy(cn_synd_mem, cn_synd_old, e->row*cir_sz, 0, cir_sz);
            // update syndrome
            vec_mod2_add(cn_flp_sel, cn_synd_old, cn_synd_new, cir_sz);
            // update syndrome memory
            vec_copy(cn_synd_new, cn_synd_mem, 0, e->row*cir_sz, cir_sz);
        }
    }

    for (int i=0; i<tm_sz*cir_sz; i++)
    {
        if (cn_synd_mem[i] != az1[i])
            int ind = i;
    }

    // G*Z1
    mod2sparse_mulvec(qc_g, enc_di_blk, gz);
    vec_clr(gz+pad_bit, cir_sz-pad_bit); // zero padding bits

    // F1*P2
    mod2sparse_mulvec(qc_f1, gz, f1z);

    // F2*P2
    mod2sparse_mulvec(qc_f2, gz, f2z);

    // A*Z1+F1*P2
    vec_mod2_add(az1, f1z, az1, tm_sz*cir_sz);

    // E*(A*Z1)
    mod2sparse_mulvec(qc_e, az1, eaz1);

    // C*Z1
    mod2sparse_mulvec(qc_c, enc_di_blk, cz1);

    // C*Z1+F2*P2
    vec_mod2_add(cz1, f2z, cz1, (bm_m-tm_sz-1)*cir_sz);

    // E*(A*Z1)+C*Z1
    vec_mod2_add(eaz1, cz1, sumz1, (bm_m-tm_sz-1)*cir_sz);

    // Z2 = F_inv*[E*(A*Z1)+C*Z1]
    mod2sparse_mulvec(qc_fi, sumz1, z2);

    // B*Z2
    mod2sparse_mulvec(qc_b, z2, bz2);

    // Z3 = BZ2 +AZ1
    vec_mod2_add(az1, bz2, z3, tm_sz*cir_sz);

    // encoded data
    vec_copy(enc_di_blk, enc_do_blk, 0, 0, hm_k);
    vec_copy(z2, enc_do_blk, 0, hm_k, (bm_m-tm_sz)*cir_sz);
    vec_copy(z3, enc_do_blk, 0, (bm_m-tm_sz-1)*cir_sz, tm_sz*cir_sz);
    vec_copy(gz, enc_do_blk, 0, (bm_n-1)*cir_sz, pad_bit);

    // removing 0 padding
    vec_copy(enc_do_blk, tx_blk, 0, 0, info_len);
    vec_copy(enc_do_blk, tx_blk, hm_k, info_len, hm_m-cir_sz+pad_bit);

    // free space
    free(az1);
    free(eaz1);
    free(cz1);
    free(sumz1);
    free(z2);
    free(bz2);
    free(z3);
    free(f1z);
    free(f2z);
    free(gz);
}

void ldpc_packet::ldpc_decoder(enum dec_model dec_mode)
{
    // add 0 padding
    vec_copy(det_blk, dec_di_blk, 0, 0, info_len);
    for (int i=0; i<pad_len; i++)
        dec_di_blk[info_len + i] = max_llr_bin;

    vec_copy(det_blk, dec_di_blk, info_len, hm_k, hm_m);
    for (int i=hm_n-1;i>=hm_n-cir_sz+pad_bit; i--)
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
            ldpc_dec_bf2(3, col_skip_itr);  // 暂时使用bf2代替bf3
        else
            ldpc_dec_bf2(3, col_skip_itr);
    }
    else if (dec_mode == LAYER)
        ldpc_dec_layer();
    else if (dec_mode == TBFDEC)
        ldpc_dec_bf2(3, col_skip_itr);  // 暂时使用bf2代替2bit_bf
    else if (dec_mode == PPBF)
    {
        // Default probability table for dv3R050N1296
        double default_p[5] = {0.0, 0.0081, 0.3, 0.7, 1.0};
        ldpc_dec_ppbf(3, default_p); // Use p_num=3 and default probability table
    }
    else if (dec_mode == PGDBF)
    {
        double p_flip = 0.5; // Single probability used by simplified PGDBF
        ldpc_dec_pgdbf(p_flip);
    }

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
    cn_synd_old = (char *)calloc(cir_sz, sizeof(*cn_synd_old));
    cn_synd_new = (char *)calloc(cir_sz, sizeof(*cn_synd_new));


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

                    vec_shift(cn_synd_sel,vn_synd_sel,cir_sz,e->shift);
                    vec_incr(vn_synd_cnt, vn_synd_sel,cir_sz);
                }


                // previous column is skipped
                if(col_skip)
                    col_skip = false;   // Column skip feature OFF
                else
                    if(col_skip_itr==0)
                        col_skip = false;// non-skip iterations
                    else
                        if((col_skip_itr>0)&&(itr<col_skip_itr))
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
                        cw_fail = 0;
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

void ldpc_packet::ldpc_dec_bf2(int p_num, int col_skip_itr)
{
    mod2entry *e;
    int synd_wt;
    int col_updt;
    int itr_updt;
    bool col_skip;

    char *dec_sb_blk; // CW soft bi
    char *cn_synd_mem;  //syndrome memory in CN order
    char *cn_synd_sel;  //selected syndrome in CN order
    char *vn_synd_sel;  //selected syndrome in VN order
    char *vn_synd_cnt;  //syndrome weight of select columns in VN order
    char *vn_hd_sel;    //current HD of selected column (from dec_do_blk) in VN order
    char *vn_sb_sel;    //soft bit of selected column (from dec_sb_blk) in VN order
    char *vn_raw_sel;    //raw data of selected column (from dec_di_blk) in VN order
    char **vn_flp_sel;  //flip flag of selected column in VN order
    int *vn_flp_col;    //column index of the vn_flp_sel
    int *vn_flp_itr;    //iteration of the vn_flp_sel
    char *cn_flp_sel;  //flip flag of selected column in CN order
    char *cn_synd_new;   // new syndrome in CN order
    char *cn_synd_old;  


    //allocate memory
    dec_sb_blk = (char *)calloc(hm_n, sizeof(*dec_sb_blk));
    cn_synd_mem = (char *)calloc(hm_m, sizeof(*cn_synd_mem));
    cn_synd_sel = (char *)calloc(cir_sz, sizeof(*cn_synd_sel));
    vn_synd_sel = (char *)calloc(cir_sz, sizeof(*vn_synd_sel));
    vn_synd_cnt = (char *)calloc(cir_sz, sizeof(*vn_synd_cnt));
    vn_hd_sel = (char *)calloc(cir_sz, sizeof(*vn_hd_sel));
    vn_sb_sel = (char *)calloc(cir_sz, sizeof(*vn_sb_sel));
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
    vec_set(dec_sb_blk, hm_n);
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

                if (i == bm_n - 1)
                {
                    for (int ii = pad_bit; ii < cir_sz; ii++)
                        vn_synd_cnt[ii] = 7;
                }

                // previous column is skipped
                if(col_skip == false)   // Column skip feature OFF
                {
                    // read raw and current HD and sb
                    vec_copy(dec_di_blk,vn_raw_sel, i*cir_sz,0,cir_sz);
                    vec_copy(dec_do_blk,vn_hd_sel, i*cir_sz,0,cir_sz);
                    vec_copy(dec_sb_blk,vn_sb_sel, i*cir_sz,0,cir_sz);

                    // pipelines
                    for (int j = p_num; j > 0; j--)
                        vec_copy(vn_flp_sel[j - 1], vn_flp_sel[j], 0, 0, cir_sz);
                    for (int j = p_num; j > 0; j--)
                    {
                        vn_flp_col[j] = vn_flp_col[j - 1];
                        vn_flp_itr[j] = vn_flp_itr[j - 1];
                    }
                    vn_flp_col[0] = i;
                    vn_flp_itr[0] = itr;

                    for (int j = 0; j< cir_sz; j++)
                    {
                        if (vn_synd_cnt[j] == 7)
                        {
                            vn_flp_sel[0][j] = 0;
                            continue;
                        }
                        if (vn_raw_sel[j] == vn_hd_sel[j]) // same
                        {
                            if (vn_sb_sel[j] == 1) // strong
                            {
                                if (vn_synd_cnt[j] >= flp_thrshd0_s[itr - 1])
                                    vn_flp_sel[0][j] = 1;
                                else
                                    vn_flp_sel[0][j] = 0;

                                if ((vn_synd_cnt[j] >= sb_thrshd0_s0[itr - 1]) && (vn_synd_cnt[j] < sb_thrshd0_s1[itr - 1]))
                                    vn_sb_sel[j] = 0;
                            }
                            else // weak
                            {
                                if (vn_synd_cnt[j] >= flp_thrshd0_w[itr - 1])
                                    vn_flp_sel[0][j] = 1;
                                else
                                    vn_flp_sel[0][j] = 0;

                                if ((vn_synd_cnt[j] < sb_thrshd0_w0[itr - 1]) || (vn_synd_cnt[j] >= sb_thrshd0_w1[itr - 1]))
                                    vn_sb_sel[j] = 1;
                            }
                        }
                        else
                        {
                            if (vn_sb_sel[j] == 1) // strong
                            {
                                if (vn_synd_cnt[j] >= flp_thrshd1_s[itr - 1])
                                    vn_flp_sel[0][j] = 1;
                                else
                                    vn_flp_sel[0][j] = 0;

                                if ((vn_synd_cnt[j] >= sb_thrshd1_s0[itr - 1]) && (vn_synd_cnt[j] < sb_thrshd1_s1[itr - 1]))
                                    vn_sb_sel[j] = 0;
                            }
                            else // weak
                            {
                                if (vn_synd_cnt[j] >= flp_thrshd1_w[itr - 1])
                                    vn_flp_sel[0][j] = 1;
                                else
                                    vn_flp_sel[0][j] = 0;

                                if ((vn_synd_cnt[j] < sb_thrshd1_w0[itr - 1]) || (vn_synd_cnt[j] >= sb_thrshd1_w1[itr - 1]))
                                    vn_sb_sel[j] = 1;
                            }
                        }
                    }

                    col_updt = vn_flp_col[p_num];
                    itr_updt = vn_flp_itr[p_num];

                    for (int j=0; j<cir_sz; j++)
                    {
                        if (vn_flp_sel[p_num][j] == 1)
                        {
                            dec_do_blk[col_updt*cir_sz+j] = (dec_do_blk[col_updt*cir_sz+j]+1)%2;
                        }
                    }

                    vec_copy(vn_sb_sel, dec_sb_blk, 0, i*cir_sz, cir_sz);
                }//non-skipped columns(flip logic)
            }// end other iteration （not 0th iteration）
            
            // update syndrome memory
            // multiple times (col_wt)
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
                if(synd_wt == 0)
                {
                    cw_fail=0;
                    cnvg_itr = itr_updt;
                    cnvg_lyr = col_updt;
                    if(fdec_early_term_en ==1)
                        break;
                }
                else
                    cw_fail = 1;
            }
        }
    }


    if((fdec_early_term_en == 0) || (cw_fail==1))
    {
        cnvg_itr = fdec_max_itr - 1;
        cnvg_lyr = bm_n - 1;
        fina_synd_wt = vec_sum(cn_synd_mem,hm_m);
    }

    free(dec_sb_blk);
    free(cn_synd_mem);
    free(cn_synd_sel);
    free(vn_synd_sel);
    free(vn_synd_cnt);
    free(vn_hd_sel);
    free(vn_sb_sel);
    free(vn_raw_sel);
    for(int i=0;i<=p_num;i++)
        free(vn_flp_sel[i]);
    free(vn_flp_sel);
    free(vn_flp_col);
    free(vn_flp_itr);
    free(cn_flp_sel);
    free(cn_synd_new);
    free(cn_synd_old);

}


void ldpc_packet::ldpc_dec_bf2_nopadding()
{
    mod2entry *e;
    int synd_wt;
    int col_updt;
    int itr_updt;
    bool col_skip;

    char *dec_sb_blk; // CW soft bi
    char *cn_synd_mem;  //syndrome memory in CN order
    char *cn_synd_sel;  //selected syndrome in CN order
    char *vn_synd_sel;  //selected syndrome in VN order
    char *vn_synd_cnt;  //syndrome weight of select columns in VN order
    char *vn_hd_sel;    //current HD of selected column (from dec_do_blk) in VN order
    char *vn_sb_sel;    //soft bit of selected column (from dec_sb_blk) in VN order
    char *vn_raw_sel;    //raw data of selected column (from dec_di_blk) in VN order
    char *vn_flp_sel;  //flip flag of selected column in VN order
    char *cn_flp_sel;  //flip flag of selected column in CN order
    char *cn_synd_new;   // new syndrome in CN order
    char *cn_synd_old;  


    //allocate memory
    dec_sb_blk = (char *)calloc(hm_n, sizeof(*dec_sb_blk));
    cn_synd_mem = (char *)calloc(hm_m, sizeof(*cn_synd_mem));
    cn_synd_sel = (char *)calloc(cir_sz, sizeof(*cn_synd_sel));
    vn_synd_sel = (char *)calloc(cir_sz, sizeof(*vn_synd_sel));
    vn_synd_cnt = (char *)calloc(cir_sz, sizeof(*vn_synd_cnt));
    vn_hd_sel = (char *)calloc(cir_sz, sizeof(*vn_hd_sel));
    vn_sb_sel = (char *)calloc(cir_sz, sizeof(*vn_sb_sel));
    vn_raw_sel = (char *)calloc(cir_sz, sizeof(*vn_raw_sel));
    vn_flp_sel = (char *)calloc(cir_sz, sizeof(vn_flp_sel));
    cn_flp_sel = (char *)calloc(cir_sz, sizeof(*cn_flp_sel));
    cn_synd_new = (char *)calloc(cir_sz, sizeof(*cn_synd_new));
    cn_synd_old = (char *)calloc(cir_sz, sizeof(*cn_synd_old));


    cw_fail = 1;
    cw_miscorr = 0;
    fina_synd_wt =0;

    vec_copy(dec_di_blk, dec_do_blk, 0,0,hm_n);
    vec_clr(cn_synd_mem,hm_m);
    vec_clr(vn_flp_sel, cir_sz);
    vec_set(dec_sb_blk, hm_n); // initial all strong

    for (int itr=0;(itr<=fdec_max_itr)&&((fdec_early_term_en==0)||(cw_fail==1));itr++)
    {
        for (int i=0; i<bm_n; i++)
        {
            if (itr==0)
            {
                vec_copy(dec_di_blk,vn_flp_sel,i*cir_sz,0,cir_sz);
                col_updt=i;
                itr_updt = itr;
            }
            else
            {
                vec_clr(vn_synd_cnt,cir_sz);
                for(e=mod2sparse_first_in_col(qc_bm,i);!mod2sparse_at_end(e);e=mod2sparse_next_in_col(e))
                {
                    vec_copy(cn_synd_mem,cn_synd_sel,e->row*cir_sz,0,cir_sz);
                    vec_shift(cn_synd_sel,vn_hd_sel,cir_sz,e->shift);
                    vec_incr(vn_synd_cnt, vn_synd_sel,cir_sz);
                }

                // previous column is skipped
                vec_copy(dec_di_blk,vn_raw_sel, i*cir_sz,0,cir_sz);
                vec_copy(dec_do_blk,vn_hd_sel, i*cir_sz,0,cir_sz);
                vec_copy(dec_sb_blk,vn_sb_sel, i*cir_sz,0,cir_sz);

                // make a flip decision
                for (int j = 0; j< cir_sz; j++)
                {
                    if (vn_raw_sel[j] == vn_hd_sel[j]) // same
                    {
                        if (vn_sb_sel[j] == 1) // strong
                        {
                            if (vn_synd_cnt[j] >= flp_thrshd0_s[itr - 1])
                                vn_flp_sel[j] = 1;
                            else
                                vn_flp_sel[j] = 0;

                            if ((vn_synd_cnt[j] >= sb_thrshd0_s0[itr - 1]) && (vn_synd_cnt[j] < sb_thrshd0_s1[itr - 1]))
                                vn_sb_sel[j] = 0;
                        }
                        else // weak
                        {
                            if (vn_synd_cnt[j] >= flp_thrshd0_w[itr - 1])
                                vn_flp_sel[j] = 1;
                            else
                                vn_flp_sel[j] = 0;

                            if ((vn_synd_cnt[j] < sb_thrshd0_w0[itr - 1]) || (vn_synd_cnt[j] >= sb_thrshd0_w1[itr - 1]))
                                vn_sb_sel[j] = 1;
                        }
                    }
                    else
                    {
                        if (vn_sb_sel[j] == 1) // strong
                        {
                            if (vn_synd_cnt[j] >= flp_thrshd1_s[itr - 1])
                                vn_flp_sel[j] = 1;
                            else
                                vn_flp_sel[j] = 0;

                            if ((vn_synd_cnt[j] >= sb_thrshd1_s0[itr - 1]) && (vn_synd_cnt[j] < sb_thrshd1_s1[itr - 1]))
                                vn_sb_sel[j] = 0;
                        }
                        else // weak
                        {
                            if (vn_synd_cnt[j] >= flp_thrshd1_w[itr - 1])
                                vn_flp_sel[j] = 1;
                            else
                                vn_flp_sel[j] = 0;

                            if ((vn_synd_cnt[j] < sb_thrshd1_w0[itr - 1]) || (vn_synd_cnt[j] >= sb_thrshd1_w1[itr - 1]))
                                vn_sb_sel[j] = 1;
                        }
                    }
                }

                itr_updt = itr;
                col_updt = i;

                for (int j=0; j<cir_sz; j++)
                {
                    if (vn_flp_sel[j] == 1)
                    {
                        dec_do_blk[col_updt*cir_sz+j] = (dec_do_blk[col_updt*cir_sz+j]+1)%2;
                    }
                }

                vec_copy(vn_sb_sel, dec_sb_blk, 0, i*cir_sz, cir_sz);
            }

            // update syndrome memory
            // multiple times (col_wt)
            for(e=mod2sparse_first_in_col(qc_bm,col_updt);!mod2sparse_at_end(e);e=mod2sparse_next_in_col(e))
            {
                vec_shift(vn_flp_sel, cn_flp_sel,cir_sz,-1*e->shift);
                vec_copy(cn_synd_mem,cn_synd_old,e->row*cir_sz,0,cir_sz);
                vec_mod2_add(cn_flp_sel, cn_synd_old, cn_synd_new, cir_sz);
                vec_copy(cn_synd_new,cn_synd_mem,0,e->row*cir_sz,cir_sz);
            }

            // H Matrix Check
            if ((itr>0) || (i==(bm_n-1)))
            {
                synd_wt = vec_sum(cn_synd_mem,hm_m);
                if(synd_wt == 0)
                {
                    cw_fail=0;
                    cnvg_itr = itr_updt;
                    cnvg_lyr = col_updt;
                    if(fdec_early_term_en ==1)
                        break;
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

    free(dec_sb_blk);
    free(cn_synd_mem);
    free(cn_synd_sel);
    free(vn_synd_sel);
    free(vn_synd_cnt);
    free(vn_hd_sel);
    free(vn_sb_sel);
    free(vn_raw_sel);
    free(vn_flp_sel);
    free(cn_flp_sel);
    free(cn_synd_new);
    free(cn_synd_old);
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
    float *cn_r_new_pre;    // New R msg in CN order of previous layer
    float *cn_app_pre;      // APP = Q + R_new in CN order of previous layer
    float *cn_app_cur;      // APP = Q + R_new in CN order of current layer
    float *cn_q_sel_cur;    // Q msg of the select circulant from current layer
    float *cn_r_old_cur;    // old R msg in CN order of previous layer
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
    dec_init = (char *)calloc(bm_n, sizeof(*dec_init));
    drop_wit_now = (char *)calloc(bm_m, sizeof(*drop_wit_now));
    vec_set(dec_init, bm_n);

    cn_c_mem = (struct cn_msg **)calloc(bm_m, sizeof(*cn_c_mem));
    for (int i = 0; i < bm_m; i++)
        cn_c_mem[i] = (struct cn_msg *)calloc(cir_sz, sizeof(*cn_c_mem[i]));
    cn_c_updt_cur = (struct cn_msg *)calloc(cir_sz, sizeof(*cn_c_updt_cur));

    cn_q_mem = (float **)calloc(bm_n, sizeof(*cn_q_mem));
    for (int i = 0; i < bm_n; i++)
        cn_q_mem[i] = (float *)calloc(cir_sz, sizeof(*cn_q_mem[i]));

    cn_r_new_pre = (float *)calloc(cir_sz, sizeof(*cn_r_new_pre));
    cn_app_pre = (float *)calloc(cir_sz, sizeof(*cn_app_pre));
    cn_app_cur = (float *)calloc(cir_sz, sizeof(*cn_app_cur));
    cn_r_old_cur = (float *)calloc(cir_sz, sizeof(*cn_r_old_cur));
    cn_q_updt_cur = (float *)calloc(cir_sz, sizeof(*cn_q_updt_cur));

    cn_q_sign = (int **)calloc(bm_n*col_wt, sizeof(*cn_q_sign));
    for (int i = 0; i<bm_n*col_wt; i++)
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
                        mem_app_start = e_pre->col / (bm_m - 1) * pad_bit % cir_sz;
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
                            cn_r_new_pre[i] = cn_c_sel_pre[i].min2_val * cn_c_sel_pre[i].sign_tot * sign_tmp;
                        else
                            cn_r_new_pre[i] = cn_c_sel_pre[i].min1_val * cn_c_sel_pre[i].sign_tot * sign_tmp;
                    }

                    // Rnew
                    // APP in CN order of previous layer
                    cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];

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
                    vn_dec_hd[i] = cn_app_pre[(i+shift_val2 + cir_sz) % cir_sz] >= 0 ? 0 : 1;
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
                            cn_r_old_cur[i] = cn_c_sel_cur[i].min2_val * cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];
                        else
                            cn_r_old_cur[i] = cn_c_sel_cur[i].min1_val * cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];
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
                        sign_tmp = (cn_q_updt_cur[i] >= 0) ? 1 : -1;
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
                    cn_c_mem[layer][i].min1_val = (float)Sat_Quan((double)cn_c_mem[layer][i].min1_val, finite_c_max, finite_c_min, finite_c_num, finite_f_num);
                    cn_c_mem[layer][i].min2_val = (float)Sat_Quan((double)cn_c_mem[layer][i].min2_val, finite_c_max, finite_c_min, finite_c_num, finite_f_num);
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
            }
            else if ((hd_updated == 0) && (layer_synd_wt == 0))
            {
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

            if ((synd_pass_cnt >= bm_m) && (hd_stable_cnt >= bm_m -1))
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

    if ((cw_fail == 1) || (ldec_early_term_en == 0))
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
    for (int i=0; i<bm_n; i++)
        free(cn_q_mem[i]);
    free(cn_q_mem);
    free(cn_r_new_pre);
    free(cn_app_pre);
    free(cn_app_cur);
    free(cn_q_sel_cur);
    free(cn_r_old_cur);
    free(cn_q_updt_cur);
    for (int i=0; i<bm_n*col_wt; i++)
        free(cn_q_sign[i]);
    free(cn_q_sign);
    free(layer_synd);
    free(cn_dec_hd);
    free(vn_dec_hd);
    free(drop_wit_now);

#ifdef _LDPC_DEBUG_DUMP
    fclose(cfp);
    fclose(sfp);
    fclose(hdfp);
    fclose(lfp);
#endif
} // ldpc_dec_layer

// gen3 layer decoder
void ldpc_packet::ldpc_dec_layer3()
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
    float *cn_r_new_pre;    // New R msg in CN order of previous layer
    float *cn_app_pre;      // APP = Q + R_new in CN order of previous layer
    float *cn_app_cur;      // APP = Q + R_new in CN order of current layer
    float *cn_q_sel_cur;    // Q msg of the select circulant from current layer
    float *cn_r_old_cur;    // old R msg in CN order of previous layer
    float *cn_q_updt_cur;   // Updated Q msg of the select circulant in current layer
    int **cn_q_sign;      // Q sign

    char *layer_synd;
    char *cn_dec_hd;
    char *vn_dec_hd;
    int hd_updated;
    int layer_synd_wt;
    int synd_pass_cnt = 0;
    int hd_stable_cnt = 0;
    
    // allocation
    dec_init = (char *)calloc(bm_n, sizeof(*dec_init));
    vec_set(dec_init, bm_n);

    cn_c_mem = (struct cn_msg **)calloc(bm_m, sizeof(*cn_c_mem));
    for (int i = 0; i<bm_m; i++)
        cn_c_mem[i] = (struct cn_msg *)calloc(cir_sz, sizeof(*cn_c_mem[i]));
    cn_c_updt_cur = (struct cn_msg *)calloc(cir_sz, sizeof(*cn_c_updt_cur));

    cn_q_mem = (float **)calloc(bm_n, sizeof(*cn_q_mem));
    for (int i = 0; i<bm_n; i++)
        cn_q_mem[i] = (float*)calloc(cir_sz, sizeof(*cn_q_mem[i]));

    cn_r_new_pre = (float *)calloc(cir_sz, sizeof(*cn_r_new_pre));
    cn_app_pre = (float *)calloc(cir_sz, sizeof(*cn_app_pre));
    cn_app_cur = (float *)calloc(cir_sz, sizeof(*cn_app_cur));
    cn_q_sel_cur = (float *)calloc(cir_sz, sizeof(*cn_q_sel_cur));
    cn_r_old_cur = (float *)calloc(cir_sz, sizeof(*cn_r_old_cur));
    cn_q_updt_cur = (float *)calloc(cir_sz, sizeof(*cn_q_updt_cur));

    cn_q_sign = (int **)calloc(bm_m*col_wt, sizeof(*cn_q_sign));
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
        for (int j=0; j<cir_sz; j++)
            cn_q_mem[i][j] = (float)llr_tbl[dec_di_blk[i*cir_sz+j]];

    // iterative decoding
    for(int itr=0;(itr<=ldec_max_itr)&&((ldec_early_term_en==0)||(cw_fail==1));itr++)
    {
        // Q sign mem index
        cir_cnt = 0;

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

                // find out the previous layer of select column
                e_pre = mod2sparse_prev_in_col(e);
                if (mod2sparse_at_end(e_pre))
                    e_pre = mod2sparse_last_in_col(qc_bm, e->col);
                
                cn_c_sel_pre = cn_c_mem[e_pre->row];    // read previous layer C msg

                // cal Rnew and APP
                for (int i = 0; i < cir_sz; i++)
                {
                    // Qmsg sign
                    sign_tmp = (cn_q_sel_pre[i] >= 0) ? 1 : -1;

                    // Rnew
                    if (cn_c_sel_pre[i].min1_pos == e->col)
                        cn_r_new_pre[i] = cn_c_sel_pre[i].min2_val * cn_c_sel_pre[i].sign_tot * sign_tmp;
                    else
                        cn_r_new_pre[i] = cn_c_sel_pre[i].min1_val * cn_c_sel_pre[i].sign_tot * sign_tmp;

                    // APP in CN order of previous layer
                    cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];

                    // Quantization
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
                    vn_dec_hd[i] = cn_app_pre[(i+shift_val2 + cir_sz) % cir_sz] >= 0 ? 0 : 1;
                }

                // CW converge check logic per circulant
                // 1. check if HD updated
                if (hd_updated == 0)
                    if (vec_cmp(dec_do_blk, vn_dec_hd, e->col*cir_sz, 0, cir_sz) == 1)
                        hd_updated = 1;

                vec_copy(vn_dec_hd, dec_do_blk, 0, e->col*cir_sz, cir_sz);

                // 2. accumulate syndrome
                vec_shift(vn_dec_hd, cn_dec_hd, cir_sz, -1 * e->shift);
                vec_mod2_add(cn_dec_hd, layer_synd, layer_synd, cir_sz);

                // calculate R_old and current Q, update current layer C and Q
                for (int i=0; i<cir_sz; i++)
                {
                    if (cn_c_sel_cur[i].min1_pos == e->col)
                        cn_r_old_cur[i] = cn_c_sel_cur[i].min2_val * cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];
                    else
                        cn_r_old_cur[i] = cn_c_sel_cur[i].min1_val * cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];

                    // Q -= Rold
                    cn_q_updt_cur[i] = cn_app_cur[i] - cn_r_old_cur[i];
                    // Quantization
                    if (finite_mode == 1)
                    {
                        cn_q_updt_cur[i] = (float)Sat_Quan((double)cn_q_updt_cur[i], finite_q_max, finite_q_min, finite_q_num, finite_f_num);
                    }

                    // update C
                    sign_tmp = (cn_q_updt_cur[i] >= 0) ? 1 : -1;
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

                    cn_q_sign[cir_cnt][i] = sign_tmp;
                }

                // update Q memory
                for (int i=0; i<cir_sz; i++)
                {
                    cn_q_mem[e->col][i] = cn_q_updt_cur[i];
                }
  
                cir_cnt++;
            } // per circulant

            // update C_MSG per layer
            for (int i=0; i<cir_sz; i++)
            {
                cn_c_mem[layer][i].min1_val = cn_c_updt_cur[i].min1_val * alpha;
                cn_c_mem[layer][i].min2_val = cn_c_updt_cur[i].min2_val * alpha;

                if (finite_mode == 1)
                {
                    cn_c_mem[layer][i].min1_val = (float)Sat_Quan((double)cn_c_mem[layer][i].min1_val, finite_c_max, finite_c_min, finite_c_num, finite_f_num);
                    cn_c_mem[layer][i].min2_val = (float)Sat_Quan((double)cn_c_mem[layer][i].min2_val, finite_c_max, finite_c_min, finite_c_num, finite_f_num);
                }

                cn_c_mem[layer][i].min1_pos = cn_c_updt_cur[i].min1_pos;
                cn_c_mem[layer][i].sign_tot = cn_c_updt_cur[i].sign_tot;
            }

            // check converage checking
            layer_synd_wt = vec_sum(layer_synd, cir_sz);
            if (hd_init == 1)
            {
                hd_stable_cnt = 0;
                synd_pass_cnt = 0;
            }
            else if ((hd_updated == 0) && (layer_synd_wt == 0))
            {
                hd_stable_cnt++;
                synd_pass_cnt++;
            }
            else
            {
                hd_stable_cnt = 0;
                synd_pass_cnt = 0;
            }

            if ((synd_pass_cnt >= bm_m) && (hd_stable_cnt >= bm_m -1))
            {
                cw_fail = 0;
                cnvg_itr = itr;
                cnvg_lyr = layer;
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
    for (int i=0; i<bm_n; i++)
        free(cn_q_mem[i]);
    free(cn_q_mem);
    free(cn_r_new_pre);
    free(cn_app_pre);
    free(cn_app_cur);
    free(cn_q_sel_cur);
    free(cn_r_old_cur);
    free(cn_q_updt_cur);
    for (int i=0; i<bm_n*col_wt; i++)
        free(cn_q_sign[i]);
    free(cn_q_sign);
    free(layer_synd);
    free(cn_dec_hd);
    free(vn_dec_hd);
} // ldpc_dec_layer3


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

void ldpc_packet::ldpc_dec_ppbf(int p_num, double *p)
{
    // Rename to PGDBF for clarity
    char *hard = dec_do_blk;  // Current hard decisions (will be modified)
    char *hard0;              // Initial hard decisions (saved)
    char *syndrome;           // Syndrome vector
    char *Ej;                  // Energy function vector
    int *flipPos;             // Positions to flip
    
    // Allocate memory
    hard0 = (char *)calloc(hm_n, sizeof(*hard0));
    syndrome = (char *)calloc(hm_m, sizeof(*syndrome));
    Ej = (char *)calloc(hm_n, sizeof(*Ej));
    flipPos = (int *)calloc(hm_n, sizeof(*flipPos)); // Maximum possible positions
    
    cw_fail = 1;
    cw_miscorr = 0;
    fina_synd_wt = 0;
    
    // Initialize: hard = initial decisions, save hard0
    vec_copy(dec_di_blk, hard, 0, 0, hm_n);    // hard = initial input
    vec_copy(dec_di_blk, hard0, 0, 0, hm_n);   // hard0 = saved initial input
    
    int iteration = 0;
    
    // Calculate initial syndrome
    mod2sparse_mulvec(qc_hm, hard, syndrome);
    int syndrome_weight = vec_sum(syndrome, hm_m);

    if (syndrome_weight == 0)
    {
        cw_fail = 0;
        cnvg_itr = 0;
        cnvg_lyr = 0;
    }
    else
    {
        while ((syndrome_weight != 0) && (iteration < fdec_max_itr))
        {
            iteration++;
        
        // Calculate energy function: Ej = mod(hard+hard0, 2) + syndrome*H
        // First part: mod(hard+hard0, 2)
        for (int j = 0; j < hm_n; j++) {
            Ej[j] = (hard[j] + hard0[j]) % 2;
        }
        
        // Second part: syndrome*H 
        // For each column j in H, calculate dot product with syndrome
        mod2entry *e;
        for (int j = 0; j < hm_n; j++) {
            int col_product = 0;
            // Iterate through non-zero entries in column j
            for (e = mod2sparse_first_in_col(qc_hm, j); 
                 !mod2sparse_at_end(e); 
                 e = mod2sparse_next_in_col(e)) {
                col_product += syndrome[e->row];
            }
            col_product = col_product % 2;
            
            // Add to energy function
            Ej[j] += col_product;
        }
        
        // Find maximum value in Ej
        int bmax = vec_max(Ej, hm_n);
        
        // Find all positions where Ej == bmax
        int num_flip_pos = vec_find(Ej, hm_n, bmax, flipPos);
        
        if (num_flip_pos > 0) {
            // Apply probabilistic flipping
            for (int k = 0; k < num_flip_pos; k++) {
                double rand_val = rand_uniform();
                if (rand_val < p[0]) { // Use first probability value
                    int pos = flipPos[k];
                    hard[pos] = (hard[pos] + 1) % 2; // Flip the bit
                }
            }
        }
        
        // Recalculate syndrome
        mod2sparse_mulvec(qc_hm, hard, syndrome);
        syndrome_weight = vec_sum(syndrome, hm_m);
        
        // Check for convergence
            if (syndrome_weight == 0)
            {
                cw_fail = 0;
                cnvg_itr = iteration - 1;
                cnvg_lyr = 0; // Not applicable for PGDBF
                break;
            }
        }

        if (cw_fail == 1)
        {
            cnvg_itr = (fdec_max_itr > 0) ? (fdec_max_itr - 1) : 0;
            cnvg_lyr = 0;
            fina_synd_wt = syndrome_weight;
        }
    }

    // Copy result to output
    vec_copy(hard, dec_do_blk, 0, 0, hm_n);
    
    // Free memory
    free(hard0);
    free(syndrome);
    free(Ej);
    free(flipPos);
    
#ifdef _LDPC_DEBUG
    if (cw_fail == 0) {
        printf("[LDPC DEBUG] PGDBF decoding success after %d iterations.\n", cnvg_itr);
    } else {
        printf("[LDPC DEBUG] PGDBF decoding failed after %d iterations, final syndrome weight: %d\n", cnvg_itr, fina_synd_wt);
    }
#endif
}

// LDPC PGDBF decoder with pipeline mechanism and padding handling
// Based on ldpc_dec_bf but implements PGDBF algorithm
void ldpc_packet::ldpc_dec_pgdbf(int p_num, int col_skip_itr, double *p_values)
{
    mod2entry *e;
    int synd_wt;
    int col_updt;
    int itr_updt;
    bool col_skip;
    char *cn_synd_mem;      // syndrome memory in CN order
    char *cn_synd_sel;      // selected syndrome in CN order
    char *vn_synd_sel;      // selected syndrome in VN order
    char *vn_energy_cnt;     // energy function values in VN order
    char *vn_hd_sel;        // current HD of selected column in VN order
    char *vn_raw_sel;       // raw data of selected column in VN order
    char **vn_flp_sel;      // flip flag of selected column in VN order (pipeline)
    int *vn_flp_col;        // column index of the vn_flp_sel
    int *vn_flp_itr;        // iteration of the vn_flp_sel
    char *cn_flp_sel;       // flip flag of selected column in CN order
    char *cn_synd_new;      // new syndrome in CN order
    char *cn_synd_old;      // old syndrome in CN order
    char *hard0_backup;     // backup of initial hard decisions
    int *flipPos;           // positions where energy is maximum
    
    // Allocate memory
    cn_synd_mem = (char *)calloc(hm_m, sizeof(*cn_synd_mem));
    cn_synd_sel = (char *)calloc(cir_sz, sizeof(*cn_synd_sel));
    vn_synd_sel = (char *)calloc(cir_sz, sizeof(*vn_synd_sel));
    vn_energy_cnt = (char *)calloc(cir_sz, sizeof(*vn_energy_cnt));
    vn_hd_sel = (char *)calloc(cir_sz, sizeof(*vn_hd_sel));
    vn_raw_sel = (char *)calloc(cir_sz, sizeof(*vn_raw_sel));
    vn_flp_sel = (char **)calloc(p_num+1, sizeof(*vn_flp_sel));
    for(int i = 0; i <= p_num; i++)
        vn_flp_sel[i] = (char *)calloc(cir_sz, sizeof(*vn_flp_sel[i]));
    vn_flp_col = (int *)calloc(p_num+1, sizeof(*vn_flp_col));
    vn_flp_itr = (int *)calloc(p_num+1, sizeof(*vn_flp_itr));
    cn_flp_sel = (char *)calloc(cir_sz, sizeof(*cn_flp_sel));
    cn_synd_new = (char *)calloc(cir_sz, sizeof(*cn_synd_new));
    cn_synd_old = (char *)calloc(cir_sz, sizeof(*cn_synd_old));
    hard0_backup = (char *)calloc(hm_n, sizeof(*hard0_backup));
    flipPos = (int *)calloc(cir_sz, sizeof(*flipPos));

    cw_fail = 1;
    cw_miscorr = 0;
    fina_synd_wt = 0;

    // Initialize: hard = initial decisions, save hard0
    vec_copy(dec_di_blk, dec_do_blk, 0, 0, hm_n);
    vec_copy(dec_di_blk, hard0_backup, 0, 0, hm_n);
    vec_clr(cn_synd_mem, hm_m);
    for(int i = 0; i < p_num; i++)
        vec_clr(vn_flp_sel[i], cir_sz);
    fdec_cyc_num = 0;
    fdec_cyc_org = 0;

    for(int itr = 0; (itr <= fdec_max_itr) && ((fdec_early_term_en == 0) || (cw_fail == 1)); itr++)
    {
        for(int i = 0; i < (itr == fdec_max_itr ? p_num : bm_n); i++)
        {
            if(itr == 0)
            {
                vec_copy(dec_di_blk, vn_flp_sel[p_num], i*cir_sz, 0, cir_sz);
                col_updt = i;
                itr_updt = itr;
                col_skip = false;
            }
            else
            {
                vec_clr(vn_energy_cnt, cir_sz);
                
                // Calculate syndrome contribution: syndrome * H
                for(e = mod2sparse_first_in_col(qc_bm, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
                {
                    int drop_start = cir_sz + 1;
                    int drop_length = 0;
                    int min_drop_ind = -1;
                    vec_copy(cn_synd_mem, cn_synd_sel, e->row*cir_sz, 0, cir_sz);

                    // Handle padding for last column
                    if(e->col == bm_n-1)
                    {
                        drop_start = (cir_sz - (e->shift) + pad_bit) % cir_sz;
                        drop_length = cir_sz - pad_bit;
                        min_drop_ind = drop_start + drop_length;
                    }
                    else
                    {
                        if(drop_col_bit_map[e->row][e->col])
                        {
                            if(e->row == bm_m - 1)
                            {
                                drop_start = pad_bit;
                                drop_length = cir_sz - pad_bit;
                            }
                            else
                            {
                                drop_start = (e->col - e->row)/(bm_m-1)*pad_bit % cir_sz;
                                drop_length = pad_bit;
                            }
                            min_drop_ind = drop_start + drop_length;
                        }
                    }

                    // Zero out padded bits
                    for(int ii = 0; ii < cir_sz; ii++)
                    {
                        if((min_drop_ind <= cir_sz && ii >= drop_start && ii < min_drop_ind)
                        || (min_drop_ind > cir_sz && (ii >= drop_start || ii < min_drop_ind % cir_sz)))
                            cn_synd_sel[ii] = 0;
                    }

                    // Apply circular shift and accumulate syndrome contribution
                    vec_shift(cn_synd_sel, vn_synd_sel, cir_sz, e->shift);
                    vec_incr(vn_energy_cnt, vn_synd_sel, cir_sz);
                }

                // Calculate energy function: Ej = mod(hard+hard0, 2) + syndrome*H
                vec_copy(dec_do_blk, vn_hd_sel, i*cir_sz, 0, cir_sz);
                vec_copy(hard0_backup, vn_raw_sel, i*cir_sz, 0, cir_sz);
                
                for(int j = 0; j < cir_sz; j++)
                {
                    // First part: mod(hard + hard0, 2)
                    int first_part = (vn_hd_sel[j] + vn_raw_sel[j]) % 2;
                    // Total energy function
                    vn_energy_cnt[j] += first_part;
                }

                // Column skip logic (same as original BF)
                if(col_skip)
                    col_skip = false;
                else
                    if(col_skip_itr == 0)
                        col_skip = false;
                    else
                        if((col_skip_itr > 0) && (itr < col_skip_itr))
                            col_skip = false;
                        else
                        {
                            col_skip = true;
                            // Make skip decision based on energy threshold
                            int max_energy = vec_max(vn_energy_cnt, cir_sz);
                            if(max_energy >= flp_thrshd1[itr-1])
                                col_skip = false;
                        }

                // PGDBF flip logic
                if(col_skip == false)
                {
                    // Pipeline management (same as original BF)
                    for(int j = p_num; j > 0; j--)
                        vec_copy(vn_flp_sel[j-1], vn_flp_sel[j], 0, 0, cir_sz);
                    for(int j = p_num; j > 0; j--)
                    {
                        vn_flp_col[j] = vn_flp_col[j-1];
                        vn_flp_itr[j] = vn_flp_itr[j-1];
                    }
                    vn_flp_col[0] = i;
                    vn_flp_itr[0] = itr;

                    // Find maximum energy and positions
                    int bmax = vec_max(vn_energy_cnt, cir_sz);
                    int num_flip_pos = vec_find(vn_energy_cnt, cir_sz, bmax, flipPos);
                    
                    // Clear flip flags
                    vec_clr(vn_flp_sel[0], cir_sz);
                    
                    // Probabilistic flipping for maximum energy positions
                    if(num_flip_pos > 0)
                    {
                        double p_flip = (itr <= fdec_max_itr && p_values != NULL) ? p_values[itr-1] : 0.5;
                        for(int k = 0; k < num_flip_pos; k++)
                        {
                            if(rand_uniform() < p_flip)
                            {
                                vn_flp_sel[0][flipPos[k]] = 1;
                            }
                        }
                    }

                    col_updt = vn_flp_col[p_num];
                    itr_updt = vn_flp_itr[p_num];

                    // Apply flips
                    for(int j = 0; j < cir_sz; j++)
                    {
                        if(vn_flp_sel[p_num][j] == 1)
                        {
                            dec_do_blk[col_updt*cir_sz+j] = (dec_do_blk[col_updt*cir_sz+j]+1)%2;
                        }
                    }
                }
            }
            fdec_cyc_org++;

            if(col_skip == false)
            {
                fdec_cyc_num++;

                // Update syndrome memory (same as original BF)
                for(e = mod2sparse_first_in_col(qc_bm, col_updt); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
                {
                    int drop_start = cir_sz + 1;
                    int drop_length = 0;
                    int min_drop_ind = -1;

                    vec_shift(vn_flp_sel[p_num], cn_flp_sel, cir_sz, -1*e->shift);
                    
                    if(col_updt == bm_n-1)
                    {
                        drop_start = (cir_sz - (e->shift) + pad_bit) % cir_sz;
                        drop_length = cir_sz - pad_bit;
                        min_drop_ind = drop_start + drop_length;
                    }
                    else
                    {
                        if(drop_col_bit_map[e->row][e->col])
                        {
                            if(e->row == bm_m-1)
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
                    
                    for(int ii = 0; ii < cir_sz; ii++)
                    {
                        if((min_drop_ind <= cir_sz && ii >= drop_start && ii < min_drop_ind)
                        || (min_drop_ind > cir_sz && (ii >= drop_start || ii < min_drop_ind % cir_sz)))
                            cn_flp_sel[ii] = 0;
                    }
                    
                    vec_copy(cn_synd_mem, cn_synd_old, e->row*cir_sz, 0, cir_sz);
                    vec_mod2_add(cn_flp_sel, cn_synd_old, cn_synd_new, cir_sz);
                    vec_copy(cn_synd_new, cn_synd_mem, 0, e->row*cir_sz, cir_sz);
                }

                // Check for convergence
                if((itr > 0) || (i == (bm_n-1)))
                {
                    synd_wt = vec_sum(cn_synd_mem, hm_m);
                    if(synd_wt == 0)
                    {
                        cw_fail = 0;
                        cnvg_itr = itr_updt;
                        cnvg_lyr = col_updt;
                        if(fdec_early_term_en == 1)
                            break;
                    }
                }
            }
        }
    }

    if((fdec_early_term_en == 0) || (cw_fail == 1))
    {
        cnvg_itr = fdec_max_itr - 1;
        cnvg_lyr = bm_n - 1;
        fina_synd_wt = vec_sum(cn_synd_mem, hm_m);
    }

    // Free memory
    free(cn_synd_mem);
    free(cn_synd_sel);
    free(vn_synd_sel);
    free(vn_energy_cnt);
    free(vn_hd_sel);
    free(vn_raw_sel);
    for(int i = 0; i <= p_num; i++)
        free(vn_flp_sel[i]);
    free(vn_flp_sel);
    free(vn_flp_col);
    free(vn_flp_itr);
    free(cn_flp_sel);
    free(cn_synd_new);
    free(cn_synd_old);
    free(hard0_backup);
    free(flipPos);

#ifdef _LDPC_DEBUG
    if (cw_fail == 0) {
        printf("[LDPC DEBUG] PGDBF decoding success after %d iterations.\n", cnvg_itr);
    } else {
        printf("[LDPC DEBUG] PGDBF decoding failed after %d iterations, final syndrome weight: %d\n", cnvg_itr, fina_synd_wt);
    }
#endif
}


void ldpc_packet::ldpc_dec_pgdbf(double p_flip)
{
    char *hard = dec_do_blk;
    char *hard0 = (char *)calloc(hm_n, sizeof(*hard0));
    char *syndrome = (char *)calloc(hm_m, sizeof(*syndrome));
    int  *energy = (int *)calloc(hm_n, sizeof(*energy));
    int  *flip_pos = (int *)calloc(hm_n, sizeof(*flip_pos));

    if (!hard0 || !syndrome || !energy || !flip_pos) {
        printf("[LDPC] Memory allocation failure in ldpc_dec_pgdbf_simple\n");
        exit(1);
    }

    cw_fail = 1;
    cw_miscorr = 0;
    fina_synd_wt = 0;
    fdec_cyc_num = 0;
    fdec_cyc_org = 0;

    vec_copy(dec_di_blk, hard, 0, 0, hm_n);
    vec_copy(dec_di_blk, hard0, 0, 0, hm_n);

    mod2sparse_mulvec(qc_hm, hard, syndrome);
    int syndrome_weight = vec_sum(syndrome, hm_m);

    int iteration = 0;

    if (syndrome_weight == 0) {
        cw_fail = 0;
        cnvg_itr = 0;
        cnvg_lyr = 0;
    } else {
        while (syndrome_weight != 0 && iteration < fdec_max_itr) {
            iteration++;
            fdec_cyc_org++;

        int max_energy = -1;
        int num_flip = 0;

        for (int j = 0; j < hm_n; j++) {
            int column_sum = 0;
            for (mod2entry *e = mod2sparse_first_in_col(qc_hm, j);
                 !mod2sparse_at_end(e);
                 e = mod2sparse_next_in_col(e)) {
                column_sum += syndrome[e->row];
            }

            int first_part = (hard[j] + hard0[j]) & 1;
            int ej = first_part + column_sum;
            energy[j] = ej;

            if (ej > max_energy) {
                max_energy = ej;
                flip_pos[0] = j;
                num_flip = 1;
            } else if (ej == max_energy) {
                flip_pos[num_flip++] = j;
            }
        }

        bool flipped = false;
        for (int i = 0; i < num_flip; i++) {
            if (rand_uniform() < p_flip) {
                int pos = flip_pos[i];
                hard[pos] ^= 1;
                flipped = true;
            }
        }
        if (flipped) fdec_cyc_num++;

            mod2sparse_mulvec(qc_hm, hard, syndrome);
            syndrome_weight = vec_sum(syndrome, hm_m);
            if (syndrome_weight == 0) {
                cw_fail = 0;
                cnvg_itr = iteration - 1;
                cnvg_lyr = 0;
                break;
            }
        }

        if (cw_fail) {
            cnvg_itr = (fdec_max_itr > 0) ? (fdec_max_itr - 1) : 0;
            cnvg_lyr = 0;
            fina_synd_wt = syndrome_weight;
        }
    }

    free(hard0);
    free(syndrome);
    free(energy);
    free(flip_pos);

#ifdef _LDPC_DEBUG
    if (!cw_fail)
        printf("[LDPC DEBUG] Simplified PGDBF decoding success after %d iterations.\n", cnvg_itr);
    else
        printf("[LDPC DEBUG] Simplified PGDBF decoding failed after %d iterations, final syndrome weight: %d\n",
               cnvg_itr, fina_synd_wt);
#endif
}

void ldpc_packet::ldpc_dec_mbf(int Fx, int threshold)
{
    char *hard = dec_do_blk;
    char *syndrome = (char *)calloc(hm_m, sizeof(*syndrome));
    mbf_candidate_t *candidates = (mbf_candidate_t *)calloc(hm_n, sizeof(*candidates));

    if (!syndrome || !candidates) {
        printf("[LDPC] Memory allocation failure in ldpc_dec_mbf\n");
        exit(1);
    }

    cw_fail = 1;
    cw_miscorr = 0;
    fina_synd_wt = 0;
    fdec_cyc_num = 0;
    fdec_cyc_org = 0;

    vec_copy(dec_di_blk, hard, 0, 0, hm_n);

    mod2sparse_mulvec(qc_hm, hard, syndrome);
    int syndrome_weight = vec_sum(syndrome, hm_m);

    int iteration = 0;

    if (syndrome_weight == 0) {
        cw_fail = 0;
        cnvg_itr = 0;
        cnvg_lyr = 0;
    } else {
        while (syndrome_weight != 0 && iteration < fdec_max_itr) {
            iteration++;
            fdec_cyc_org++;

        int candidate_count = 0;

        for (int j = 0; j < hm_n; j++) {
            int energy_val = 0;
            for (mod2entry *e = mod2sparse_first_in_col(qc_hm, j);
                 !mod2sparse_at_end(e);
                 e = mod2sparse_next_in_col(e)) {
                energy_val += ((int)syndrome[e->row] * 2) - 1;
            }

            if (energy_val > threshold) {
                candidates[candidate_count].idx = j;
                candidates[candidate_count].energy = energy_val;
                candidate_count++;
            }
        }

        bool flipped = false;
        if (candidate_count > 0) {
            int max_flips = (Fx > 0) ? Fx : candidate_count;
            int flips_to_apply = candidate_count;

            if (candidate_count > max_flips) {
                qsort(candidates, candidate_count, sizeof(*candidates), compare_mbf_candidate_desc);
                flips_to_apply = max_flips;
            }

            for (int i = 0; i < flips_to_apply; i++) {
                int pos = candidates[i].idx;
                hard[pos] ^= 1;
                flipped = true;
            }
        }

        if (flipped)
            fdec_cyc_num++;

            mod2sparse_mulvec(qc_hm, hard, syndrome);
            syndrome_weight = vec_sum(syndrome, hm_m);
            if (syndrome_weight == 0) {
                cw_fail = 0;
                cnvg_itr = iteration - 1;
                cnvg_lyr = 0;
                break;
            }
        }

        if (cw_fail) {
            cnvg_itr = (fdec_max_itr > 0) ? (fdec_max_itr - 1) : 0;
            cnvg_lyr = 0;
            fina_synd_wt = syndrome_weight;
        }
    }

    free(syndrome);
    free(candidates);

#ifdef _LDPC_DEBUG
    if (!cw_fail)
        printf("[LDPC DEBUG] MBF decoding success after %d iterations.\n", cnvg_itr);
    else
        printf("[LDPC DEBUG] MBF decoding failed after %d iterations, final syndrome weight: %d\n",
               cnvg_itr, fina_synd_wt);
#endif
}

/*
ldpc_packet::ldpc_dec_bf_ibex(s_ldpc_decoder_input ldpc_decoder_input, s_ldpc_decoder_parameters ldpc_decoder_parameters, s_h_matrix h_matrix)
{
    int VERBOSITY = 0;
    int MAX_ERROR_COUNT = 4095;
    int i;
    int j;
    int k;
    int m;
    int iteration = 0;
    int clock_cycles = 0;
    int syndrome_weight;
    int weight;
    int decode_mode;
    int decode_mode2;
    bool post_process = 0;
    bool post_process_2 = 0;
    bool do_post_flipped = 0;
    bool do_post_unflipped = 0;
    bool give_up = 0;
    bool finished = 0;
    bool post_processing = 0;
    bool qc_parity_en;
    bool disable_update;
    bool hamming_weight_le_circ_thr;
    bool hamming_weight_lt_circ_thr;
    bool post_trigger;
    bool post_trigger2;
    bool flipped_prev;
    bool flipped;
    bool be_aggressive;
    bool prng_post_process;
    bool prng_post_process2;
    s_hard_codeword hard_codeword;
    s_hard_codeword soft_codeword;
    s_variable_nodes vn;
    s_check_node cn;
    s_check_node cn_shifted;
    s_likelihood_levels likelihood_levels;
    s_256_bits prng_256;
    s_256_bits prng_512;
    int prng_init[32];
    int syndrome_weight_delayed;

    int syndrome_weight_r;
    bool do_not_use_this_bit;
    bool look;

    prng_init[31] = 0x083d;
    prng_init[30] = 0x3214;
    prng_init[29] = 0xa8a1;
    prng_init[28] = 0x5327;
    prng_init[27] = 0x71bc;
    prng_init[26] = 0x3edb;
    prng_init[25] = 0xba50;
    prng_init[24] = 0xc946;
    prng_init[23] = 0x4c9a;
    prng_init[22] = 0x0b73;
    prng_init[21] = 0xef18;
    prng_init[20] = 0x31ee;
    prng_init[19] = 0xf2b2;
    prng_init[18] = 0xd98f;
    prng_init[17] = 0x89f9;
    prng_init[16] = 0x9375;
    prng_init[15] = 0x043d;
    prng_init[14] = 0x6214;
    prng_init[13] = 0xa8a1;
    prng_init[12] = 0x5d27;
    prng_init[11] = 0x71cc;
    prng_init[10] = 0x2edb;
    prng_init[9]  = 0xb550;
    prng_init[8]  = 0xc646;
    prng_init[7]  = 0x8c9a;
    prng_init[6]  = 0x1b73;
    prng_init[5]  = 0xef08;
    prng_init[4]  = 0x30ee;
    prng_init[3]  = 0xf7b2;
    prng_init[2]  = 0xda8f;
    prng_init[1]  = 0x49f9;
    prng_init[0]  = 0x9365;

    for (j=0; j< h_matrix.cols; j++)
    {
        for (k=0; k<h_matrix.bits; k++)
        {
            hard_codeword.c[j].b[k] = ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard;
            vn.c[j].b[k].bit_hard = hard_codeword.c[j].b[k];
            soft_codeword.c[j].b[k] = ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable;
            vn.c[j].b[k].flipped = 0;
        }
    }

    CN = f_check_nodes(h_matrix, hard_codeword);
    if (VERBOSITY > 0)
    {
        printf("[LDPC DEBUG] Starting BF decoding with max %d iterations.\n", fdec_max_itr);
        f_print_hard_codeword(hard_codeword, h_matrix.cols, h_matrix.bits);
        printf("### DECODER C++: H MATRIX:\n");
        f_print_h_matrix(h_matrix);
        printf("### DECODER C++: AFTER SYNDROME CHECK, CODEWORD HAS CHECK NODES:\n");
        f_print_check_nodes(cn, h_matrix.rows, h_matrix.bits);
    }

    for (i=0; i<h_matrix.rows; i++)
        for (j=0; j<h_matrix.bits; j++)
            cn_shifted.r[i].b[k] = 0;
    syndrome_weight = f_check_node_weight(h_matrix, cn);
    clock_cycles = (2 * h_matrix.rows) + 1;
    if (syndrome_weight == 0)
        finished = 1;
    ldpc_decoder_output.syndrome_weight_before = syndrome_weight;
    syndrome_weight_delayed = syndrome_weight;
    ldpc_decoder_output.early_termination = 0;

    if (ldpc_decoder_parameters.early_terminate_dis == 0)
    {
        int early_term_thr;
        if (ldpc_decoder_input.nand_strobes == 0)
        {
            if      (h_matrix.rows <=6)
                early_term_thr = ldpc_decoder_parameters.early_terminate_thr[0][0];
            else if (h_matrix.rows == 7)
                early_term_thr = ldpc_decoder_parameters.early_terminate_thr[0][1];
            else if (h_matrix.rows == 8)
                early_term_thr = ldpc_decoder_parameters.early_terminate_thr[0][2];
            else if (h_matrix.rows == 9)
                early_term_thr = ldpc_decoder_parameters.early_terminate_thr[0][3];
            else if (h_matrix.rows == 10)
                early_term_thr = ldpc_decoder_parameters.early_terminate_thr[0][4];
            else if (h_matrix.rows == 11)
                early_term_thr = ldpc_decoder_parameters.early_terminate_thr[0][5];
            else
                early_term_thr = ldpc_decoder_parameters.early_terminate_thr[0][6];
        }    
        else
        {
            if      (h_matrix.rows <=6)
                early_term_thr = ldpc_decoder_parameters.early_terminate_thr[1][0];
            else if (h_matrix.rows == 7)
                early_term_thr = ldpc_decoder_parameters.early_terminate_thr[1][1];
            else if (h_matrix.rows == 8)
                early_term_thr = ldpc_decoder_parameters.early_terminate_thr[1][2];
            else if (h_matrix.rows == 9)
                early_term_thr = ldpc_decoder_parameters.early_terminate_thr[1][3];
            else if (h_matrix.rows == 10)
                early_term_thr = ldpc_decoder_parameters.early_terminate_thr[1][4];
            else if (h_matrix.rows == 11)
                early_term_thr = ldpc_decoder_parameters.early_terminate_thr[1][5];
            else
                early_term_thr = ldpc_decoder_parameters.early_terminate_thr[1][6];
        }

        if (syndrome_weight >= early_term_thr)
        {
            finished = 1;
            ldpc_decoder_output.early_termination = 1;
            if (VERBOSITY > 0 && syndrome_weight >= early_term_thr)
                printf("### LDPC DECODER C++ EARLY TERMINATION: %5d\n", syndrome_weight);
        }
    }
    
}
*/