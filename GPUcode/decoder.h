#pragma once
#include <chrono>
#include <cstdint>
#include "codec_support.h"
#include "codeword.h"
#include "h_matrix.h"

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
    void print_stats(void);
};

class decoder : public codec_support {
  public:
    static uint16_t prng_init[32];

    struct ldpc_decoder_params {
        static const uint8_t VN_BITS = 8;
        bool post_process_en = true;
        uint16_t syndrome_weight_thr_qc = 16;
        uint16_t syndrome_weight_thr_post = 48;
        uint16_t likelihood_thr = 128;
        uint16_t likelihood_init_coef_all[8][4] = {
            {0, 7, 6, 5}, {0, 7, 5, 5}, {0, 6, 5, 4}, {0, 5, 5, 4}, {0, 4, 4, 3}, {0, 3, 3, 2}, {0, 4, 3, 2}, {0, 4, 3, 2},
        };
        uint16_t likelihood_init_fraction[3] = {0, 6, 6}; // strobe 5: {0,6,6}, strobe 7: {0,8,4}
        uint16_t post_ratio = 12;
    };

    ldpc_decoder_params ldpc_decoder_parameters;
    const h_matrix &h_matrix_ref;
    decoder(h_matrix &h_matrix_ref, int nand_strobes = 5, int post_ratio = 12);
    ldpc_decoder_output decode_planar(const ldpc_decoder_input &input);

  private:
    void update_vn_planar(uint64_t *vn_flipped, const uint8_t *weight_word, unsigned char *vn_likelihood, const s_likelihood_levels &likelihood_levels, const bool *prng_512, bool post_trigger1, bool post_trigger2, int bit_offset);
    s_likelihood_levels compute_likelihood_levels(int strobes, decoder::ldpc_decoder_params ldpc_decoder_parameters, int syndrome_weight, int rows);
    s_512_bits lfsr_512_bit(s_512_bits data_in);
};

struct logger {
    static const uint64_t AVG_WINDOW_SIZE = 1;
    static const uint64_t MAX_ITER = 2048;
    static uint64_t ITER_LIMIT;
    static void log_elapsed_time(const std::chrono::steady_clock::time_point &start_time, const std::chrono::steady_clock::time_point &end_time);
    static void print_accumulated_stats(uint64_t accumulated_cw_count, ldpc_decoder_output &decode_stats);
    static void print_accumulated_stats(uint64_t accumulated_cw_count, decoder_output_acc &decode_stats);
    static void print_iter_stats(uint32_t *iter_info);
};
