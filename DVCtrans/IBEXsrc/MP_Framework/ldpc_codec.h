// LDPC packet inherit from CH packet
// 0. LDPC configuration/clean
// 1. LDPC packet allocation/clean
// 2. LDPC encoder
// 3. LDPC bit-flipping decoder
// 4. LDPC min-sum decoder

#ifndef _LDPC_PACKET_H
#define _LDPC_PACKET_H

#include "mod2sparse.h"
#include "transceiver.h"

enum
{
    LDPC_MAX_ROWS = 17,
    LDPC_MAX_COLS = 84,
    LDPC_MAX_PAYLOAD_COLS = 67,
    LDPC_MAX_CIRC_BITS = 512,
    LDPC_PMS_LUT_SIZE = 6,
};

enum dec_model
{
    SKIP = 0,
    BF_P0,
    BF_P3 = 3,
    BF_G2 = 4,
    BF_IBEX = 5,
    BF_IBEX_RTL_CN = 6,
    LAYER = 8,
    LAYER_G2 = 9,
};

struct cn_msg
{
    float min1_val;
    float min2_val;
    int min1_pos;
    char sign_tot;
};

struct s_h_matrix
{
    int rows;
    int cols;
    int bits;
    int column_weight;
    int bytes_of_userdata;
    int bytes_of_parity;
    int extra_bytes_of_parity;
    int extra_bytes_of_userdata;
    int unused_bytes_of_parity;
    int unused_bytes_of_userdata;
    int extra_bits_of_parity;
    int extra_bits_of_userdata;
    int min_rows;
    int max_rows;
    int delta[LDPC_MAX_ROWS];
    int element[LDPC_MAX_ROWS][LDPC_MAX_COLS];
    int first_element[LDPC_MAX_ROWS];
    int last_element[LDPC_MAX_ROWS];
    int wraparound[LDPC_MAX_ROWS];
    int wrap_base[LDPC_MAX_ROWS];
    int wrap_num_deltas[LDPC_MAX_ROWS];
    int row_weight[LDPC_MAX_ROWS];
    int col_weight[LDPC_MAX_COLS];
    bool occupied[LDPC_MAX_ROWS][LDPC_MAX_COLS];
    bool fade[LDPC_MAX_ROWS][LDPC_MAX_COLS];
    bool mask[LDPC_MAX_COLS][LDPC_MAX_CIRC_BITS];
    bool parity_column[LDPC_MAX_COLS];
    int bits_in_last_column;
};

struct s_hard_codeword_column
{
    bool b[512];
};

struct s_hard_codeword
{
    s_hard_codeword_column c[LDPC_MAX_COLS];
    unsigned int errors_at_level_and_weight[4][5];
    unsigned int correct_at_level_and_weight[4][5];
    double probability_of_error_at_level_and_weight[4][5];
};

struct s_codeword_bit
{
    bool bit_hard;
    bool bit_questionable;
    bool bit_questionable2;
    short int level;
    short int syndrome_weight;
    bool bit_corrected;
    bool bit_is_error;
};

struct s_codeword_column
{
    s_codeword_bit b[512];
};

struct s_codeword
{
    s_codeword_column c[LDPC_MAX_COLS];
};

struct s_ldpc_decoder_parameters
{
    int post_process_en;
    int syndrome_weight_thr_qc;
    int syndrome_weight_thr_post;
    int likelihood_thr;
    int likelihood_init_coef_all[8][4];
    int likelihood_init_coef[4];
    int likelihood_init_fraction[3];
    int soft_bit_table[8];
    int post_ratio;
    int likelihood_map[8];
    bool early_terminate_dis;
    int early_terminate_thr[2][8];
    bool questionable_sense;
};

struct s_ldpc_decoder_input
{
    int nand_strobes;
    int soft_bits;
    int post_iteration;
    int iteration_limit;
    bool syndrome_cal_only;
    int errors_in_userdata;
    int errors;
    float rber;
    bool verbose;
    s_codeword corrupted_codeword;
};

struct s_ldpc_decoder_output
{
    int failure;
    int errors_in_userdata;
    int errors_in_codeword;
    int iterations;
    int clock_cycles;
    int syndrome_weight_before;
    int syndrome_weight_after;
    int early_termination;
    int col_cnt; // 0-based column index where decoding converged; -1 if not converged or no column processed
    s_hard_codeword corrected_codeword;
};

struct s_variable_node_bit
{
    short int likelihood;
    bool bit_hard;
    bool bit_questionable;
    bool bit_questionable2;
    short int level;
    bool flipped;
};

struct s_variable_node_column
{
    s_variable_node_bit b[512];
};

struct s_variable_nodes
{
    s_variable_node_column c[LDPC_MAX_COLS];
};

struct s_check_node_row
{
    bool b[512];
};

struct s_check_nodes
{
    s_check_node_row r[LDPC_MAX_ROWS];
};

struct s_256_bits
{
    bool b[256];
};

struct s_512_bits
{
    bool b[512];
};

struct s_likelihood_levels
{
    int min;
    int max;
    int flip_thr;
    int weak;
    int strong;
    int level[4];
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

    // layer decoder scheduler constraints
    int rdec_cmem_cont_thrshd;
    int rdec_hdmem_cont_thrshd;

    // QC code
    mod2sparse *qc_bm;  // base matrix
    mod2sparse *qc_hm;  // QC H matrix
    mod2sparse *qc_a, *qc_b, *qc_c, *qc_d, *qc_e, *qc_fi; // encoder matrix

    // IBEX
    s_h_matrix h_matrix;
    s_ldpc_decoder_input ldpc_decoder_input;
    s_ldpc_decoder_output ldpc_decoder_output;
    s_ldpc_decoder_parameters ldpc_decoder_parameters;
    int total_cir;

    // decoder config
    int fdec_max_itr;
    int fdec_early_term_en;
    int ldec_max_itr;
    int ldec_early_term_en;
    int *flp_thrshd0;
    int *flp_thrshd1;
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

    float alpha;
    float alpha_pms[LDPC_PMS_LUT_SIZE];
    float beta_pms[LDPC_PMS_LUT_SIZE];
    float point1;
    float point2;
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

    // SDLite LLR override (DV config passthrough)
    int reg_sdlite_llr_config;
    int reg_sdlite_llr0;
    int reg_sdlite_llr1;
    int reg_sdlite_llr2;
    int reg_sdlite_llr3;

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

    void ldpc_config(int, int, int, int, int);
    void ldpc_dec_config(int, int, int, float, const float *, const float *,
                         float, float, int, int, int, int, int, int, int, int,
                         int, int);
    void ldpc_gen_gm();
    void ldpc_clean();
    void ldpc_ibex_phck(s_h_matrix);

    // IBEX core
    void ldpc_ibex_encoder();
    void ldpc_ibex_input(int, int, int, int);
    void ldpc_ibex_parameters(int,
                              int,
                              int,
                              int,
                              int,
                              int,
                              int,
                              int,
                              int,
                              unsigned int,
                              unsigned int,
                              unsigned int,
                              unsigned int,
                              unsigned int,
                              unsigned int,
                              unsigned int);
    s_hard_codeword f_ldpc_encode(s_hard_codeword, s_h_matrix);
    void f_print_h_matrix(s_h_matrix);
    void f_print_hard_codeword(s_hard_codeword, int, int);
    s_check_nodes f_check_nodes(s_h_matrix, s_hard_codeword);
    int f_check_node_weight(s_h_matrix, s_check_nodes);
    s_likelihood_levels f_likelihood_levels(int, s_ldpc_decoder_parameters, int, int);
    s_256_bits f_256_bit_lfsr(s_256_bits);
    s_512_bits f_512_bit_lfsr(s_512_bits);
    void f_print_s_256_bits(s_256_bits);
    void f_print_s_512_bits(s_512_bits);
    void f_print_check_nodes(s_check_nodes, int, int);
    void f_print_check_nodes_shifted(s_check_nodes, s_h_matrix, int);
    int f_update_vn_post(int, int, int, int, bool, bool, bool, int, bool);
    
    // print H matrix
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

    // LDPC IBEX BF decoder
    void ldpc_dec_bf_ibex(s_ldpc_decoder_input, s_ldpc_decoder_parameters, s_h_matrix);
    void ldpc_dec_bf_ibex_rtl_cn(s_ldpc_decoder_input, s_ldpc_decoder_parameters, s_h_matrix);
    // void ldpc_dec_bf_ibex_2bit(s_ldpc_decoder_input, s_ldpc_decoder_parameters, s_h_matrix);

    // LDPC layer decoder
    void ldpc_dec_layer();
    void ldpc_dec_layer2();

    // LDPC PMS decoder
    void ldpc_dec_pms();
    int ldpc_pms_ind(float, float);
    // Skip LDPC decoder
    void ldpc_dec_skip();

};

#endif // _LDPC_PACKET_H
