#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "finite_lib.h"
#include "ldpc_codec.h"
#include "mod2convert.h"
#include "mod2dense.h"
#include "mod2sparse.h"
#include "vec_op.h"

void ldpc_packet::ldpc_rd_phck(char *pchk_file)
{
    mod2entry *e;
    FILE *fp;
    int col_shift;
    qc_bm = mod2sparse_allocate(bm_m, bm_n);
    qc_hm = mod2sparse_allocate(hm_m, hm_n);

    printf("[LDPC] Porting H matrix from %s ...\n", pchk_file);

    fp = fopen(pchk_file, "r");

    if (fp == NULL)
    {
        printf("[LDPC] Error opening file %s \n", pchk_file);
        exit(0);
    }

    for (int i = 0; i < bm_m; i++)
    {
        for (int j = 0; j < bm_n; j++)
        {
            fscanf(fp, "%d", &col_shift);

            // check submatrix T (identity matrix)
            if ((i < tm_sz) && (j > (bm_n - tm_sz - 1)))
            {
                if (((i != (j + tm_sz - bm_n)) && (col_shift >= 0)) || ((i == (j + tm_sz - bm_n)) && (col_shift != 0)))
                {
                    printf("[LDPC] Error: Submatrix T is not identity matrix!\n");
                    exit(0);
                }
            }

            // insert mod2entry to base matrix
            if (col_shift >= 0)
            {
                e = mod2sparse_insert(qc_bm, i, j);
                e->shift = col_shift;

                for (int k=0; k < cir_sz; k++)
                    mod2sparse_insert(qc_hm, i*cir_sz + k, j*cir_sz + (k + col_shift) % cir_sz);
            }
        }
    }

    fclose(fp);
    printf("[LDPC] H matrix porting ready!\n");
}

// Generate G matrices from H
void ldpc_packet::ldpc_gen_gm()
{
    FILE *fp;
    mod2sparse *qc_ac, *qc_bd, *qc_te;
    mod2sparse *qc_exb, *qc_f;
    mod2dense *qc_f_d, *qc_fi_d;
    int *qc_a_cols, *qc_b_cols, *qc_t_cols;
    int *qc_a_rows, *qc_c_rows;      

    qc_ac = mod2sparse_allocate(hm_m, hm_k);
    qc_bd = mod2sparse_allocate(hm_m, (bm_m-tm_sz)*cir_sz);
    qc_te = mod2sparse_allocate(hm_m, tm_sz*cir_sz);

    qc_a_cols = (int *)calloc((bm_n-bm_m)*cir_sz, sizeof(*qc_a_cols));
    qc_b_cols = (int *)calloc((bm_m-tm_sz)*cir_sz, sizeof(*qc_b_cols));
    qc_t_cols = (int *)calloc(tm_sz*cir_sz, sizeof(*qc_t_cols));
    qc_a_rows = (int *)calloc(tm_sz*cir_sz, sizeof(*qc_a_rows));
    qc_c_rows = (int *)calloc((bm_m-tm_sz)*cir_sz, sizeof(*qc_c_rows));

    printf("[LDPC] Generating A/B/C/D/E matrices from H ...\n");

    for (int i=0; i< (bm_n - bm_m) *cir_sz; i++)
        qc_a_cols[i] = i;
    for (int i=0; i<(bm_m-tm_sz)*cir_sz; i++)
        qc_b_cols[i] = i +  (bm_n - bm_m) *cir_sz;
    for (int i=0; i<tm_sz*cir_sz; i++)
        qc_t_cols[i] = i + (bm_n-tm_sz)*cir_sz;

    for (int i=0; i<tm_sz*cir_sz; i++)
        qc_a_rows[i] = i;
    for (int i=0; i<(bm_m-tm_sz)*cir_sz; i++)
        qc_c_rows[i] = i + tm_sz*cir_sz;

    // split column first
    mod2sparse_copycols(qc_hm, qc_ac, qc_a_cols);
    mod2sparse_copycols(qc_hm, qc_bd, qc_b_cols);
    mod2sparse_copycols(qc_hm, qc_te, qc_t_cols);

    // then split row
    mod2sparse_copyrows(qc_ac, qc_a, qc_a_rows);
    mod2sparse_copyrows(qc_ac, qc_c, qc_c_rows);
    mod2sparse_copyrows(qc_bd, qc_b, qc_a_rows);
    mod2sparse_copyrows(qc_bd, qc_d, qc_c_rows);
    mod2sparse_copyrows(qc_te, qc_e, qc_c_rows);
    printf("[LDPC] A/B/C/D/E matrices ready!\n");

#ifdef _LDPC_DUMP
    char MH[50] = "H_matrix.txt";
    char MA[50] = "A_matrix.txt";
    char MB[50] = "B_matrix.txt";
    char MC[50] = "C_matrix.txt";
    char MD[50] = "D_matrix.txt";
    char ME[50] = "E_matrix.txt";

    printf("[LDPC] Dump H matrix to file %s\n", MH);
    fp = fopen(MH, "w");
    mod2sparse_print(fp, qc_hm);
    fclose(fp);
    printf("[LDPC] Dump A matrix to file %s\n", MA);
    fp = fopen(MA, "w");
    mod2sparse_print(fp, qc_a);
    fclose(fp);
    printf("[LDPC] Dump B matrix to file %s\n", MB);
    fp = fopen(MB, "w");
    mod2sparse_print(fp, qc_b);
    fclose(fp);
    printf("[LDPC] Dump C matrix to file %s\n", MC);
    fp = fopen(MC, "w");
    mod2sparse_print(fp, qc_c);
    fclose(fp);
    printf("[LDPC] Dump D matrix to file %s\n", MD);
    fp = fopen(MD, "w");
    mod2sparse_print(fp, qc_d);
    fclose(fp);
    printf("[LDPC] Dump E matrix to file %s\n", ME);
    fp = fopen(ME, "w");
    mod2sparse_print(fp, qc_e);
    fclose(fp);
#endif    

    // generate inverse F matrix (F=E*B+D)
    printf("[LDPC] Generating inverse F matrix from H matrix ...\n");

    qc_exb = mod2sparse_allocate((bm_m-tm_sz)*cir_sz, (bm_m-tm_sz)*cir_sz);
    qc_f = mod2sparse_allocate((bm_m-tm_sz)*cir_sz, (bm_m-tm_sz)*cir_sz);

    mod2sparse_multiply(qc_e, qc_b, qc_exb);
    mod2sparse_add(qc_exb, qc_d, qc_f);

    qc_f_d = mod2dense_allocate((bm_m-tm_sz)*cir_sz, (bm_m-tm_sz)*cir_sz);
    qc_fi_d = mod2dense_allocate((bm_m-tm_sz)*cir_sz, (bm_m-tm_sz)*cir_sz);

    mod2sparse_to_dense(qc_f, qc_f_d);
    mod2dense_invert(qc_f_d, qc_fi_d);
    mod2dense_to_sparse(qc_fi_d, qc_fi);

    printf("[LDPC] Inverse F matrix generated!\n");

#ifdef _LDPC_FI_DUMP
    char MFi[50] = "./output/Fi_matrix.txt";
    printf("[LDPC] Dump inverse F matrix to file %s\n", MFi);
    fp = fopen(MFi, "w");
    mod2sparse_print(fp, qc_fi);
    fclose(fp);
#endif

    // free up all internals
    mod2sparse_free(qc_ac);
    mod2sparse_free(qc_bd);
    mod2sparse_free(qc_te);
    mod2sparse_free(qc_exb);
    mod2sparse_free(qc_f);
    mod2dense_free(qc_f_d);
    mod2dense_free(qc_fi_d);
    free(qc_a_cols);
    free(qc_b_cols);
    free(qc_t_cols);
    free(qc_a_rows);
    free(qc_c_rows);
} // ldpc_gen_gm

// LDPC code configuration
void ldpc_packet::ldpc_config(int m, int n, int sc, int st, int wt, char *pchk_file)
{
    // QC matrix sizes
    bm_m = m;
    bm_n = n;
    bm_k = n - m;
    tm_sz = st;
    cir_sz = sc;

    col_wt = wt;
    hm_m = m * sc;
    hm_n = n * sc;
    hm_k = hm_n - hm_m;
    pad_len = hm_k - info_len;

    printf("[LDPC] Configuring LDPC code with m=%d, n=%d, sc=%d, st=%d, col_wt=%d\n", 
           bm_m, bm_n, cir_sz, tm_sz, col_wt);

    // read parity check matrix
    ldpc_rd_phck(pchk_file);

    // G matrices
    qc_a = mod2sparse_allocate(tm_sz * cir_sz, (bm_n - bm_m) * cir_sz);
    qc_b = mod2sparse_allocate(tm_sz * cir_sz, (bm_m - tm_sz) * cir_sz);
    qc_c = mod2sparse_allocate((bm_m - tm_sz) * cir_sz, (bm_n - bm_m) * cir_sz);
    qc_d = mod2sparse_allocate((bm_m - tm_sz) * cir_sz, (bm_m - tm_sz) * cir_sz);
    qc_e = mod2sparse_allocate((bm_m - tm_sz) * cir_sz, tm_sz * cir_sz);
    qc_fi = mod2sparse_allocate((bm_m - tm_sz) * cir_sz, (bm_m - tm_sz) * cir_sz);

    // generate G matrices
    ldpc_gen_gm();

#ifdef _LDPC_DUMP
        rdec_cmem_cont_thrshd = 10;
        rdec_hdmem_cont_thrshd = 6;

        // print_hm();
#endif
} // ldpc_config

// LDPC decoder config
void ldpc_packet::ldpc_dec_config(int max_fdec_itr, int fdec_col_skip, int max_ldec_itr, float dec_alpha, int fin_mode, int fin_q_num, int fin_r_num, int fin_f_num)
{
    // syndrome weight
    init_synd_wt_min = hm_m;
    init_synd_wt_max = 0;

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
}

void ldpc_packet::ldpc_pckt_alloc()
{
    ch_pckt_alloc();

    usr_blk = (char *)calloc(info_len, sizeof(*usr_blk));
    enc_di_blk = (char *)calloc(hm_k, sizeof(*enc_di_blk));
    enc_do_blk = (char *)calloc(hm_n, sizeof(*enc_do_blk));
    dec_di_blk = (char *)calloc(hm_n, sizeof(*dec_di_blk));
    dec_do_blk = (char *)calloc(hm_n, sizeof(*dec_do_blk));
    dec_blk = (char *)calloc(blk_len, sizeof(*dec_blk));

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
    char *az1, *eaz1, *cz1, *sumz1, *z2, *bz2, *z3;

    az1 = (char *)calloc(tm_sz*cir_sz, sizeof(*az1));
    eaz1 = (char *)calloc((bm_m-tm_sz)*cir_sz, sizeof(*eaz1));
    cz1 = (char *)calloc((bm_m-tm_sz)*cir_sz, sizeof(*cz1));
    sumz1 = (char *)calloc((bm_m-tm_sz)*cir_sz, sizeof(*sumz1));
    z2 = (char *)calloc((bm_m-tm_sz)*cir_sz, sizeof(*z2));
    bz2 = (char *)calloc(tm_sz*cir_sz, sizeof(*bz2));
    z3 = (char *)calloc(tm_sz*cir_sz, sizeof(*z3));

    // padding 0s
    vec_copy(usr_blk, enc_di_blk, 0, 0, info_len);
    for (int i=0; i<pad_len; i++)
        enc_di_blk[hm_k-pad_len+i] = 0;

    // A*Z1
    mod2sparse_mulvec(qc_a, enc_di_blk, az1);
    // E*(A*Z1)
    mod2sparse_mulvec(qc_e, az1, eaz1);
    // C*Z1
    mod2sparse_mulvec(qc_c, enc_di_blk, cz1);
    // E*(A*Z1)+C*Z1
    vec_mod2_add(eaz1, cz1, sumz1, (bm_m-tm_sz)*cir_sz);
    // Z2 = F_inv*[E*(A*Z1)+C*Z1]
    mod2sparse_mulvec(qc_fi, sumz1, z2);
    // B*Z2
    mod2sparse_mulvec(qc_b, z2, bz2);
    // Z3 = BZ2 +AZ1
    vec_mod2_add(az1, bz2, z3, tm_sz*cir_sz);

    // encoded data
    vec_copy(enc_di_blk, enc_do_blk, 0, 0, hm_k);
    vec_copy(z2, enc_do_blk, 0, hm_k, (bm_m-tm_sz)*cir_sz);
    vec_copy(z3, enc_do_blk, 0, (bm_n-tm_sz)*cir_sz, tm_sz*cir_sz);

    // removing 0 padding
    vec_copy(enc_do_blk, tx_blk, 0, 0, info_len);
    vec_copy(enc_do_blk, tx_blk, hm_k, info_len, hm_m);

    // free space
    free(az1);
    free(eaz1);
    free(cz1);
    free(sumz1);
    free(z2);
    free(bz2);
    free(z3);
}

void ldpc_packet::ldpc_decoder(enum dec_model dec_mode)
{
    // add 0 padding
    vec_copy(det_blk, dec_di_blk, 0, 0, info_len);
    for (int i=0; i < pad_len; i++)
        dec_di_blk[info_len + i] = max_llr_bin;

    vec_copy(det_blk, dec_di_blk, info_len, hm_k, hm_m);

    if (dec_mode == SKIP)
        ldpc_dec_skip();
    else if (dec_mode == BF_P0)
        ldpc_dec_bf(0, col_skip_itr);
    else if (dec_mode == BF_P3)
        ldpc_dec_bf(3, col_skip_itr);
    else if (dec_mode == BF_G2)
        ldpc_dec_bf2();
    else if (dec_mode == LAYER)
        ldpc_dec_layer();

    //remove padding
    vec_copy(dec_do_blk, dec_blk, 0, 0, info_len);
    vec_copy(dec_do_blk, dec_blk, hm_k, info_len, hm_m);

    // check error bit number
    dec_err_num = 0;
    cor_err_num = 0;
    for (int i=0; i<blk_len; i++)
    {
        if (dec_blk[i] != tx_blk[i])
            dec_err_num++;
        if (dec_blk[i] != det_blk[i])
            cor_err_num++;
    }

    // mis-correction case
    if ((cw_fail == 0) && (dec_err_num != 0))
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


    for (int itr = 0;(itr <= fdec_max_itr) && ((fdec_early_term_en == 0)||(cw_fail == 1)); itr++)
    {
        for (int i = 0; i < (itr == fdec_max_itr ? p_num : bm_n); i++)
        {
            if (itr == 0)
            {
                vec_copy(dec_di_blk,vn_flp_sel[p_num],i*cir_sz,0,cir_sz);
                col_updt=i;
                itr_updt = itr;
                col_skip=false;
            }
            else
            {
                vec_clr(vn_synd_cnt,cir_sz);
                for(e=mod2sparse_first_in_col(qc_bm,i); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
                {
                    // read syndrome
                    vec_copy(cn_synd_mem,cn_synd_sel,e->row*cir_sz,0,cir_sz);
                    // barrel shift
                    vec_shift(cn_synd_sel, vn_synd_sel, cir_sz, e->shift);
                    // increment vn_synd_cnt
                    vec_incr(vn_synd_cnt, vn_synd_sel, cir_sz);
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
                            for(int j = 0; j < cir_sz; j++)
                            {
                                if(vn_synd_cnt[j] >= flp_thrshd1[itr - 1])
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
                    // barrel shift
                    vec_shift(vn_flp_sel[p_num], cn_flp_sel,cir_sz,-1*e->shift);
                    // read old syndrome
                    vec_copy(cn_synd_mem,cn_synd_old,e->row*cir_sz,0,cir_sz);
                    // update new syndrome
                    vec_mod2_add(cn_flp_sel, cn_synd_old, cn_synd_new, cir_sz);
                    // update syndrome memory
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
}

void ldpc_packet::ldpc_dec_bf2()
{
    mod2entry *e;
    int synd_wt;
    int col_updt;
    int itr_updt;

    char *dec_sb_blk; // CW soft bi
    char *cn_synd_mem;  //syndrome memory in CN order
    char *cn_synd_sel;  //selected syndrome in CN order
    char *vn_synd_sel;  //selected syndrome in VN order
    char *vn_synd_cnt;  //syndrome weight of select columns in VN order
    char *vn_sb_sel;    //soft bit of selected column (from dec_sb_blk) in VN order
    char *vn_hd_sel;    //current HD of selected column (from dec_do_blk) in VN order
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
    vn_sb_sel = (char *)calloc(cir_sz, sizeof(*vn_sb_sel));
    vn_hd_sel = (char *)calloc(cir_sz, sizeof(*vn_hd_sel));
    vn_raw_sel = (char *)calloc(cir_sz, sizeof(*vn_raw_sel));
    vn_flp_sel = (char *)calloc(cir_sz, sizeof(*vn_flp_sel));
    cn_flp_sel = (char *)calloc(cir_sz, sizeof(*cn_flp_sel));
    cn_synd_old = (char *)calloc(cir_sz, sizeof(*cn_synd_old));
    cn_synd_new = (char *)calloc(cir_sz, sizeof(*cn_synd_new));

    // Initialize decoder
    cw_fail = 1;
    cw_miscorr = 0;
    fina_synd_wt =0;
    vec_copy(dec_di_blk, dec_do_blk, 0,0,hm_n);
    vec_clr(cn_synd_mem,hm_m);
    vec_clr(vn_flp_sel, cir_sz);
    vec_set(dec_sb_blk, hm_n);

    for(int itr=0;(itr<=fdec_max_itr)&&((fdec_early_term_en==0)||(cw_fail==1));itr++)
    {
        for(int i = 0;i < bm_n; i++)
        {
            if(itr==0)
            {
                vec_copy(dec_di_blk, vn_flp_sel, i*cir_sz, 0, cir_sz);
                col_updt=i;
                itr_updt = itr;
            }
            else
            {
                vec_clr(vn_synd_cnt,cir_sz);
                for(e=mod2sparse_first_in_col(qc_bm,i);!mod2sparse_at_end(e);e=mod2sparse_next_in_col(e))
                {
                    vec_copy(cn_synd_mem, cn_synd_sel, e->row*cir_sz, 0, cir_sz);
                    vec_shift(cn_synd_sel, vn_synd_sel, cir_sz, e->shift);
                    vec_incr(vn_synd_cnt, vn_synd_sel, cir_sz);
                }
                
                // read raw and current HD and sb
                vec_copy(dec_di_blk, vn_raw_sel, i*cir_sz,0,cir_sz);
                vec_copy(dec_do_blk, vn_hd_sel, i*cir_sz,0,cir_sz);
                vec_copy(dec_sb_blk, vn_sb_sel, i*cir_sz,0,cir_sz);

                for (int j = 0; j < cir_sz; j++)
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
                        dec_do_blk[col_updt * cir_sz + j] = (dec_do_blk[col_updt * cir_sz + j] + 1) % 2;
                    }
                }

                vec_copy(vn_sb_sel, dec_sb_blk, 0, i*cir_sz, cir_sz);
            }
            
            // update syndrome memory
            // multiple times (col_wt)
            for(e=mod2sparse_first_in_col(qc_bm,col_updt);!mod2sparse_at_end(e);e=mod2sparse_next_in_col(e))
            {
                // barrel shift
                vec_shift(vn_flp_sel, cn_flp_sel, cir_sz, -1 * e->shift);
                // read old syndrome
                vec_copy(cn_synd_mem, cn_synd_old, e->row * cir_sz, 0, cir_sz);
                // update new syndrome
                vec_mod2_add(cn_flp_sel, cn_synd_old, cn_synd_new, cir_sz);
                // update syndrome memory
                vec_copy(cn_synd_new, cn_synd_mem, 0, e->row * cir_sz, cir_sz);
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
            }
        }
    }

    if((fdec_early_term_en == 0) || (cw_fail==1))
    {
        cnvg_itr = fdec_max_itr - 1;
        cnvg_lyr = bm_n - 1;
        fina_synd_wt = vec_sum(cn_synd_mem, hm_m);
    }

    free(dec_sb_blk);
    free(cn_synd_mem);
    free(cn_synd_sel);
    free(vn_synd_sel);
    free(vn_synd_cnt);
    free(vn_sb_sel);
    free(vn_hd_sel);
    free(vn_raw_sel);
    free(vn_flp_sel);
    free(cn_flp_sel);
    free(cn_synd_new);
    free(cn_synd_old);
} // ldpc_dec_bf2

// gen3 layer decoder
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
