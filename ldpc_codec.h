// LDPC packet inherit from CH packet
// 0. LDPC configuration/clean
// 1. LDPC packet allocation/clean
// 2. LDPC encoder
// 3. LDPC bit-flipping decoder
// 4. LDPC min-sum decoder

#ifndef _LDPC_PACKET_H
#define _LDPC_PACKET_H

enum dec_model
{
    SKIP = 0,
    BF_P0,
    BF_P3 = 3,
    BF_G2 = 4,
    LAYER = 8,
    TBFDEC = 9,
    PPBF = 10,
};

struct cn_msg
{
    float min1_val;
    float min2_val;
    int min1_pos;
    char sign_tot;
};

struct ldpc_packet : ch_packet
{
    // QC matrix config
    int cir_sz; // size of the circulant matrix
    int bm_m; // number of rows in the base matrix
    int bm_n; // number of columns in the base matrix
    int bm_k;
    int tm_sz; // size of the submatrix T (indentity matrix)
    int col_wt; // column weight

    int hm_m; // number of rows in the H matrix
    int hm_n; // number of columns in the H matrix
    int hm_k; // user length of the H matrix
    
    int pad_len; // length of the padding bits
    int pad_bit; // padding bit num

    // layer decoder scheduler constraints
    int rdec_cmem_cont_thrshd;
    int rdec_hdmem_cont_thrshd;

    // QC code
    mod2sparse *qc_bm;  // base matrix
    mod2sparse *qc_hm;  // QC H matrix
    mod2sparse *qc_a, *qc_b, *qc_c, *qc_d, *qc_e, *qc_fi, *qc_g, *qc_f1, *qc_f2; // encoder matrix

    // Drop col info
    int **drop_col;
    int **drop_col_bit_map;
    int *wit_drop;
    int k_val;
    int awon;

    // decoder config
    int fdec_max_itr;
    int tbfdec_max_itr;
    int fdec_early_term_en;
    int ldec_max_itr;
    int ldec_early_term_en;
    int *flp_thrshd0;
    int *flp_thrshd1;
    int *tbbf_thrshd0;
    int *tbbf_thrshd1;
    int *flp_thrshd0_s;
    int *flp_thrshd0_w;
    int *flp_thrshd1_s;
    int *flp_thrshd1_w;
    int *sb_thrshd0_s0;
    int *sb_thrshd0_s1;
    int *sb_thrshd0_w0;
    int *sb_thrshd0_w1;
    int *sb_thrshd1_s0;
    int *sb_thrshd1_s1;
    int *sb_thrshd1_w0;
    int *sb_thrshd1_w1;  

    int *fpd_flp_thrshd0_s;
    int *fpd_flp_thrshd0_w;
    int *fpd_flp_thrshd1_s;
    int *fpd_flp_thrshd1_w;
    int *fpd_sb_thrshd0_s0;
    int *fpd_sb_thrshd0_s1;
    int *fpd_sb_thrshd0_w0;
    int *fpd_sb_thrshd0_w1;
    int *fpd_sb_thrshd1_s0;
    int *fpd_sb_thrshd1_s1;
    int *fpd_sb_thrshd1_w0;
    int *fpd_sb_thrshd1_w1;      

    float alpha;
    int finite_mode;
    int finite_q_num;
    int finite_r_num;
    int finite_c_num;
    int finite_f_num;
    float finite_q_max;
    float finite_q_min;
    float finite_r_max;
    float finite_r_min;
    float finite_c_max;
    float finite_c_min;
    int reg_fp_flg;

    // data block
    char *usr_blk;
    char *enc_di_blk;
    char *enc_do_blk;
    char *dec_di_blk;
    char *dec_do_blk;
    char *dec_blk;

    // CW status
    int cw_fail;
    int cw_miscorr;
    int cor_err_num;
    int dec_err_num;
    int col_skip_itr;
    int init_synd_wt;
    int init_synd_wt_max;
    int init_synd_wt_min;
    int fina_synd_wt;
    int cnvg_itr;
    int cnvg_lyr;
    int fdec_cyc_num;
    int fdec_cyc_org;

    // alloc/cleanup packet
    void ldpc_pckt_alloc();
    void ldpc_pckt_clean();

    // QC-LDPC config & clean up
    void ldpc_config(int, int, int, int, int, int, int, char *, char *);
    void ldpc_dec_config(int, int, int, int, float, int, int, int, int, int);
    void ldpc_rd_phck(char *, char *);
    void ldpc_gen_gm();
    void ldpc_clean();

    void print_hm();
    //syndrome check
    int ldpc_synd(char *);
    // LDPC enc
    void ldpc_encoder();
    // LDPC dec
    void ldpc_decoder(enum dec_model);
    // LDPC BF decoder
    void ldpc_dec_bf(int p_num, int col_skip_itr);
    // LDPC BF Gen2
    void ldpc_dec_bf2(int p_num, int col_skip_itr);
    // LDPC layer decoder
    void ldpc_dec_layer();
    // LDPC 2bit BF decoder
    void ldpc_dec_2bit_bf(int p_num, int col_skip_itr);
    // LDPC enhance BF decoder with flip memory
    void ldpc_dec_bf3(int p_num, int col_skip_itr);
    // LDPC PPBF decoder (Probabilistic Parallel Bit-Flipping)
    void ldpc_dec_ppbf(int max_ite);
    // Skip LDPC decoder
    void ldpc_dec_skip();

};

#endif // _LDPC_PACKET_H
