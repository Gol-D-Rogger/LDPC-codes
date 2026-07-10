#pragma once

#include <cuda_fp16.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <ctime>

#include "decoder.h"

#define USE_CN_MSG

struct row_info {
    uint16_t actual_idx;
    uint16_t mask_weight;
    uint16_t mask_offset;
    uint16_t element;
};

struct device_global_info_struct {
    short cols;
    short rows;
    short active_rows_in_col[LDPC_N];
    row_info elem_info_trnsp[LDPC_N][5];
    char e_pre[LDPC_M][LDPC_N];

    uint16_t nand_strobes;
    uint16_t soft_bits;
    uint16_t post_iteration;
    uint16_t iteration_limit;

    char VN_BITS;
    bool post_process_en;
    // uint16_t syndrome_weight_thr_qc;
    // uint16_t syndrome_weight_thr_post;
    // uint16_t likelihood_thr;
    // uint16_t post_ratio;
    // uint16_t likelihood_init_coef_all[8][4];
    // uint16_t likelihood_init_fraction[3];

    float error_region_prob[7];
    float vref[7];
    uint8_t max, flip_thr;

    bool prng_verilog_mode;

    uint16_t bits;
    int bytes_of_userdata;
    int bytes_of_parity;
    uint16_t element[LDPC_M][LDPC_N];
    char fade[LDPC_M][LDPC_N];
    char occupied[LDPC_M][LDPC_N];
    uint16_t mask[80][32];
    char extra_bytes_of_parity;
    int extra_bits_of_parity;
    int extra_bytes_of_userdata;
    int finite_mode;
    int finite_q_num;
    int finite_c_num;
    int finite_f_num;
    float finite_q_max;
    float finite_q_min;
    float finite_r_max;
    float finite_r_min;
    float finite_c_max;
    float finite_c_min;
    float alpha;
    float llr_table[8];
    float awgn_sigma;
};

struct check_nodes_gpu {
    static constexpr int WORD_COUNT = 8;
    static constexpr int WORD_SIZE = 64;
    static constexpr int MAX_ROW_COUNT = LDPC_M;

    uint64_t rows[MAX_ROW_COUNT][WORD_COUNT];
};

struct cn_msg {
    __half min1_val;
    __half min2_val;
    char min1_pos;
    char sign_tot;
};

struct variable_nodes_gpu {
    static constexpr int WORD_COUNT = 8;
    static constexpr int WORD_SIZE = 64;
    static constexpr int MAX_COL_COUNT = LDPC_N;
    static constexpr int MAX_ROW_COUNT = LDPC_M;

#ifdef USE_CN_MSG
    cn_msg cn_c_mem[MAX_COL_COUNT][WORD_COUNT * WORD_SIZE];
#else
    __half cn_c_mem_min1_val[MAX_COL_COUNT][WORD_COUNT * WORD_SIZE];
    __half cn_c_mem_min2_val[MAX_COL_COUNT][WORD_COUNT * WORD_SIZE];
    uint8_t cn_c_mem_min1_pos[MAX_COL_COUNT][WORD_COUNT * WORD_SIZE];
#endif
    __half cn_q_mem[MAX_COL_COUNT][WORD_COUNT * WORD_SIZE];
    codeword dec_do_blk;
};

struct decoder_workspace {
    static constexpr int MAX_COL_COUNT = LDPC_N;

    variable_nodes_gpu vn;
    uint32_t cn_q_sign[5 * MAX_COL_COUNT][16];
    cn_msg cn_c_updt_cur[512];
};

struct OptimizedSharedMemory {
    static constexpr int MAX_COL_COUNT = LDPC_N;

    uint16_t syndrome_weight;
    bool finished;

    uint8_t min;
    uint8_t likelihood_levels[4];

    row_info row_elem_info[5];
    __half cn_app_pre[512];
    uint16_t layer_synd[32], vn_dec_hd[32];
    int dec_init[MAX_COL_COUNT];
};

void copy_decoder_info_to_device(const device_global_info_struct &host_decoder_info);
int execute_gpu_kernel(const uint64_t TOTAL_CODEWORDS, const uint64_t MAX_FAILURE_COUNT, decoder_input_cw *test_cw = NULL);
