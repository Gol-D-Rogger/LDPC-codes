#pragma once

#include <chrono>
#include <cstdint>
#include <iomanip>
#include <sstream>

#include "codeword.h"
#include "h_matrix.h"
#include "codec_support.h"

struct ldpc_decoder_input {
    uint64_t nand_strobes = 0;
    uint64_t soft_bits = 0;
    uint64_t post_iteration = 1024;
    uint64_t iteration_limit = 2048;
    decoder_input_cw corrupted_codeword;

    ldpc_decoder_input();
    ldpc_decoder_input(uint64_t nand_strobes, uint64_t soft_bits, uint64_t iteration_limit, uint64_t post_iteration);
    void clear_cw(void);
};

struct cn_msg_cpu {
    float min1_val;
    float min2_val;
    int min1_pos;
    char sign_tot;
};

struct ldpc_decoder_output {
    uint32_t failure = 0;
    uint32_t total_errors = 0;
    uint32_t iterations = 0;
    uint32_t clock_cycles = 0;
    uint32_t syndrome_weight_before = 0;
    uint32_t syndrome_weight_after = 0;
    uint32_t net_acc = 0;

    void print_stats(void);
};

struct decoder_output_acc {
    uint64_t failure = 0;
    uint64_t total_errors = 0;
    uint64_t iterations = 0;
    uint64_t clock_cycles = 0;
    uint64_t syndrome_weight_before = 0;
    uint64_t syndrome_weight_after = 0;
    uint64_t net_acc = 0;
};

class decoder : public codec_support {
  public:
    static uint16_t prng_init[32];

    struct ldpc_decoder_params {
        static const uint8_t VN_BITS = 8;
        bool post_process_en = 1;
        int finite_mode = 1;
        int finite_q_num = 9;
        int finite_f_num = 3;
        int finite_c_num = 7;
        int llr_tot_bit = 7;
        int llr_frac_bit = 3;
        float finite_q_max = 31.875;
        float finite_q_min = -31.875;
        float finite_r_max = 7.875;
        float finite_r_min = -7.875;
        float finite_c_max = 7.875;
        float finite_c_min = 0;
        float alpha = 0.625;
        float llr_table[8] = {7.875000, 6.625000, 3.750000, 1.250000, -7.875000, -6.625000, -3.750000, -1.250000};
    };

    ldpc_decoder_params ldpc_decoder_parameters;
    const h_matrix &h_matrix_ref;

    decoder(h_matrix &h_matrix_ref, int nand_strobes = 5, int post_ratio = 12);
    ldpc_decoder_output decode_planar(const ldpc_decoder_input &input);

  private:
    void update_vn_planar(uint64_t *vn_flipped, const uint8_t *weight_word, unsigned char *vn_likelihood, const s_likelihood_levels &likelihood_levels, const bool *prng_512, const bool post_trigger1, const bool post_trigger2, const int bit_offset);
    s_likelihood_levels compute_likelihood_levels(int strobes, decoder::ldpc_decoder_params ldpc_decoder_parameters, int syndrome_weight, int rows);
    s_512_bits lfsr_512_bit(s_512_bits data_in);
    void update_node(float *cn_q_sel_pre, cn_msg_cpu *cn_c_sel_pre, float *cn_r_new_pre, float *cn_app_pre, int layer_pre, int col);
    void update_hd(char *dec_init, int col, int layer, int layer_pre, float *cn_app_pre, float *cn_app_cur, char *vn_dec_hd, int &hd_updated, char *dec_do_blk, char *cn_dec_hd, char *layer_synd);
    void update_node_next(float *cn_r_old_cur, cn_msg_cpu *cn_c_sel_cur, float *cn_q_updt_cur, float *cn_app_pre, float *cn_app_cur, int finite_mode, int cir_cnt, int layer, int col, int **cn_q_sign, float **cn_q_mem, cn_msg_cpu *cn_c_updt_cur);
};

struct logger {
    static const auto AVG_WINDOW_SIZE = 1;
    static const auto MAX_ITER = 2048;
    static uint64_t ITER_LIMIT;

    static void log_elapsed_time(
        const std::chrono::steady_clock::time_point &start_time, 
        const std::chrono::steady_clock::time_point &end_time
    );
    
    static void print_accumulated_stats(uint64_t accumulated_cw_count, ldpc_decoder_output &decode_stats);
    static void print_accumulated_stats(uint64_t accumulated_cw_count, decoder_output_acc &decode_stats);
    static void print_iter_stats(uint32_t *iter_info);
};
