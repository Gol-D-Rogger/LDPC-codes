#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "decoder.h"
#include "finite_lib.h"
#include "vec_op.h"

#define COMPARE_GPU

uint16_t decoder::prng_init[] = {0x9365, 0x49f9, 0xda8f, 0xf7b2, 0x30ee, 0xef08, 0x1b73, 0x8c9a, 0xc646, 0xb550, 0x2edb, 0x71cc, 0x5d27, 0xa8a1, 0x6214, 0x043d, 0x9375, 0x89f9, 0xd98f, 0xf2b2, 0x31ee, 0xef18, 0x0b73, 0x4c9a, 0xc946, 0xba50, 0x3edb, 0x71bc, 0x5327, 0xa8a1, 0x3214, 0x083d};

uint64_t logger::ITER_LIMIT = logger::MAX_ITER;

ldpc_decoder_input::ldpc_decoder_input() {
    return;
}

ldpc_decoder_input::ldpc_decoder_input(uint64_t nand_strobes, uint64_t soft_bits, uint64_t iteration_limit, uint64_t post_iteration) {
    this->nand_strobes = nand_strobes;
    this->soft_bits = soft_bits;
    this->iteration_limit = iteration_limit;
    this->post_iteration = post_iteration;
}

void ldpc_decoder_input::clear_cw(void) {
    corrupted_codeword.clear();
}

void ldpc_decoder_output::print_stats(void) {
    std::cout << "===== Codeword Decoding Stats =====\n";
    std::cout << "failure: " << failure << "\n";
    std::cout << "iterations: " << iterations << "\n";
    std::cout << "clock_cycles: " << clock_cycles << "\n";
    std::cout << "syndrome_weight_before: " << syndrome_weight_before << "\n";
    std::cout << "syndrome_weight_after: " << syndrome_weight_after << "\n";
}

inline void print_512float(float *a) {
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 16; j++)
            printf("%f ", a[i * 16 + j]);
        printf("\n");
    }
}

inline void print_512bit(char *a) {
    for (int i = 0; i < 16; i++) {
        uint32_t tmp = 0;
        for (int j = 0; j < 32; j++)
            tmp |= (a[i * 32 + j] << j);
        printf("%08x ", tmp);
    }
    printf("\n");
}

inline void print_512cn_msg(cn_msg_cpu *a) {
    printf("print min1_val\n");
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 16; j++)
            printf("%f ", a[i * 16 + j].min1_val);
        printf("\n");
    }

    printf("print min2_val\n");
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 16; j++)
            printf("%f ", a[i * 16 + j].min2_val);
        printf("\n");
    }

    printf("print min1_pos\n");
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 16; j++)
            printf("%2d ", a[i * 16 + j].min1_pos);
        printf("\n");
    }

    printf("print tot_sign\n");
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 16; j++)
            printf("%2d ", a[i * 16 + j].sign_tot);
        printf("\n");
    }
}

decoder::decoder(h_matrix &h_matrix_ref, int nand_strobes, int post_ratio) : h_matrix_ref(h_matrix_ref) {
    return;
}

void decoder::update_node(float *cn_q_sel_pre, cn_msg_cpu *cn_c_sel_pre, float *cn_r_new_pre, float *cn_app_pre, int layer_pre, int col) {
    for (int i = 0; i < h_matrix_ref.bits; i++) {
        const int sign_tmp = (cn_q_sel_pre[i] >= 0) ? 1 : -1;

        if (cn_c_sel_pre[i].min1_pos == col)
            cn_r_new_pre[i] = cn_c_sel_pre[i].min2_val * cn_c_sel_pre[i].sign_tot * sign_tmp;
        else
            cn_r_new_pre[i] = cn_c_sel_pre[i].min1_val * cn_c_sel_pre[i].sign_tot * sign_tmp;

        if (h_matrix_ref.extra_bytes_of_parity == 0) {
            cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
        } else if (h_matrix_ref.occupied[layer_pre][col] && (layer_pre < (h_matrix_ref.rows - 1))) {
            cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
        } else if (h_matrix_ref.occupied[layer_pre][col] && (layer_pre == (h_matrix_ref.rows - 1)) && h_matrix_ref.mask[col][(i + h_matrix_ref.element[layer_pre][col]) % h_matrix_ref.bits]) {
            cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
        } else if (h_matrix_ref.fade[layer_pre][col] && !h_matrix_ref.mask[col][(i + h_matrix_ref.element[layer_pre][col]) % h_matrix_ref.bits]) {
            cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
        } else {
            cn_app_pre[i] = cn_q_sel_pre[i];
        }

        if (ldpc_decoder_parameters.finite_mode == 1)
            cn_app_pre[i] = (float)Sat_Quan((double)cn_app_pre[i], ldpc_decoder_parameters.finite_q_max, ldpc_decoder_parameters.finite_q_min, ldpc_decoder_parameters.finite_q_num, ldpc_decoder_parameters.finite_f_num);
    }
}

void decoder::update_hd(char *dec_init, int col, int layer, int layer_pre, float *cn_app_pre, float *cn_app_cur, 
                        char *vn_dec_hd, int &hd_updated, char *dec_do_blk, char *cn_dec_hd, char *layer_synd) {
    int shift_val1;
    int shift_val2;

    if (dec_init[col] == 1) {
        shift_val1 = h_matrix_ref.element[layer][col];
        shift_val2 = 0;
        dec_init[col] = 0;
    } else {
        shift_val1 = -1 * h_matrix_ref.element[layer_pre][col] + h_matrix_ref.element[layer][col];
        shift_val2 = -1 * h_matrix_ref.element[layer_pre][col];
    }

    for (int i = 0; i < h_matrix_ref.bits; i++) {
        cn_app_cur[i] = cn_app_pre[(i + shift_val1 + h_matrix_ref.bits) % h_matrix_ref.bits];
        vn_dec_hd[i] = (cn_app_pre[(i + shift_val2 + h_matrix_ref.bits) % h_matrix_ref.bits] >= 0) ? 0 : 1;
    }

    if (hd_updated == 0) {
        if (vec_cmp(dec_do_blk, vn_dec_hd, col * h_matrix_ref.bits, 0, h_matrix_ref.bits) == 1)
            hd_updated = 1;
    }

    vec_copy(vn_dec_hd, dec_do_blk, 0, col * h_matrix_ref.bits, h_matrix_ref.bits);
    vec_shift(vn_dec_hd, cn_dec_hd, h_matrix_ref.bits, -1 * h_matrix_ref.element[layer][col]);

    if (h_matrix_ref.extra_bytes_of_parity == 0) {
        vec_mod2_add(cn_dec_hd, layer_synd, layer_synd, h_matrix_ref.bits);
    } else {
        if (h_matrix_ref.occupied[layer][col] && (layer == (h_matrix_ref.rows - 1))) {
            for (int i = 0; i < h_matrix_ref.bits; i++) {
                if (!h_matrix_ref.mask[col][(i + h_matrix_ref.element[layer][col]) % h_matrix_ref.bits])
                    cn_dec_hd[i] = 0;
            }
        }

        if (h_matrix_ref.fade[layer][col]) {
            for (int i = 0; i < h_matrix_ref.bits; i++) {
                if (h_matrix_ref.mask[col][(i + h_matrix_ref.element[layer][col]) % h_matrix_ref.bits])
                    cn_dec_hd[i] = 0;
            }
        }

        vec_mod2_add(cn_dec_hd, layer_synd, layer_synd, h_matrix_ref.bits);
    }
}

void decoder::update_node_next(float *cn_r_old_cur, cn_msg_cpu *cn_c_sel_cur, float *cn_q_updt_cur, float *cn_app_pre, float *cn_app_cur, int finite_mode, int cir_cnt, int layer, int col, int **cn_q_sign, float **cn_q_mem, cn_msg_cpu *cn_c_updt_cur) {
    (void)cn_app_pre;

    for (int i = 0; i < h_matrix_ref.bits; i++) {
        if (cn_c_sel_cur[i].min1_pos == col)
            cn_r_old_cur[i] = cn_c_sel_cur[i].min2_val * cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];
        else
            cn_r_old_cur[i] = cn_c_sel_cur[i].min1_val * cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];

        if (h_matrix_ref.extra_bytes_of_parity == 0) {
            cn_q_updt_cur[i] = cn_app_cur[i] - cn_r_old_cur[i];
        } else {
            if (h_matrix_ref.occupied[layer][col] && (layer == (h_matrix_ref.rows - 1)) && (!h_matrix_ref.mask[col][(i + h_matrix_ref.element[layer][col]) % h_matrix_ref.bits]))
                cn_r_old_cur[i] = 0;

            if (h_matrix_ref.fade[layer][col] && h_matrix_ref.mask[col][(i + h_matrix_ref.element[layer][col]) % h_matrix_ref.bits])
                cn_r_old_cur[i] = 0;

            cn_q_updt_cur[i] = cn_app_cur[i] - cn_r_old_cur[i];
        }

        if (finite_mode == 1)
            cn_q_updt_cur[i] = (float)Sat_Quan((double)cn_q_updt_cur[i], ldpc_decoder_parameters.finite_q_max, ldpc_decoder_parameters.finite_q_min, ldpc_decoder_parameters.finite_q_num, ldpc_decoder_parameters.finite_f_num);

        int sign_tmp;
        float val_tmp;
        if (h_matrix_ref.extra_bits_of_parity == 0) {
            sign_tmp = (cn_q_updt_cur[i] >= 0) ? 1 : -1;
            val_tmp = cn_q_updt_cur[i] * sign_tmp;
        } else if (h_matrix_ref.occupied[layer][col] && (layer == (h_matrix_ref.rows - 1)) && (!h_matrix_ref.mask[col][(i + h_matrix_ref.element[layer][col]) % h_matrix_ref.bits])) {
            sign_tmp = 1;
            val_tmp = 100000;
        } else if (h_matrix_ref.fade[layer][col] && h_matrix_ref.mask[col][(i + h_matrix_ref.element[layer][col]) % h_matrix_ref.bits]) {
            sign_tmp = 1;
            val_tmp = 100000;
        } else {
            sign_tmp = (cn_q_updt_cur[i] >= 0) ? 1 : -1;
            val_tmp = cn_q_updt_cur[i] * sign_tmp;
        }

        cn_c_updt_cur[i].sign_tot *= sign_tmp;
        if (val_tmp < cn_c_updt_cur[i].min1_val) {
            cn_c_updt_cur[i].min2_val = cn_c_updt_cur[i].min1_val;
            cn_c_updt_cur[i].min1_val = val_tmp;
            cn_c_updt_cur[i].min1_pos = col;
        } else if (val_tmp < cn_c_updt_cur[i].min2_val) {
            cn_c_updt_cur[i].min2_val = val_tmp;
        }

        cn_q_sign[cir_cnt][i] = sign_tmp;
    }

    for (int i = 0; i < h_matrix_ref.bits; i++) {
        cn_q_mem[col][i] = cn_q_updt_cur[i];
    }
}

ldpc_decoder_output decoder::decode_planar(const ldpc_decoder_input &decoder_input) {
    int col;
    int layer_pre;
    int cir_cnt;
    int hd_init;
    int hd_updated;
    int layer_synd_wt = 0;
    int synd_pass_cnt = 0;
    int hd_stable_cnt = 0;
    int cw_fail = 1;
    int ldec_early_term_en = 1;
    int clock_cycles = 0;

    char *dec_init = (char *)calloc(h_matrix_ref.cols, sizeof(*dec_init));
    char *layer_synd = (char *)calloc(h_matrix_ref.bits, sizeof(*layer_synd));
    char *cn_dec_hd = (char *)calloc(h_matrix_ref.bits, sizeof(*cn_dec_hd));
    char *vn_dec_hd = (char *)calloc(h_matrix_ref.bits, sizeof(*vn_dec_hd));
    char *dec_di_blk = (char *)calloc(h_matrix_ref.bits * h_matrix_ref.cols, sizeof(*dec_di_blk));
    char *dec_do_blk = (char *)calloc(h_matrix_ref.bits * h_matrix_ref.cols, sizeof(*dec_do_blk));
    vec_set(dec_init, h_matrix_ref.cols);

    cn_msg_cpu **cn_c_mem = (cn_msg_cpu **)calloc(h_matrix_ref.rows, sizeof(*cn_c_mem));
    for (int i = 0; i < h_matrix_ref.rows; i++)
        cn_c_mem[i] = (cn_msg_cpu *)calloc(h_matrix_ref.bits, sizeof(*cn_c_mem[i]));

    cn_msg_cpu *cn_c_updt_cur = (cn_msg_cpu *)calloc(h_matrix_ref.bits, sizeof(*cn_c_updt_cur));
    cn_msg_cpu *cn_c_sel_cur = NULL;
    cn_msg_cpu *cn_c_sel_pre = NULL;
    float **cn_q_mem = (float **)calloc(h_matrix_ref.cols, sizeof(*cn_q_mem));
    for (int i = 0; i < h_matrix_ref.cols; i++)
        cn_q_mem[i] = (float *)calloc(h_matrix_ref.bits, sizeof(*cn_q_mem[i]));

    float *cn_q_sel_pre = NULL;
    float *cn_r_new_pre = (float *)calloc(h_matrix_ref.bits, sizeof(*cn_r_new_pre));
    float *cn_app_pre = (float *)calloc(h_matrix_ref.bits, sizeof(*cn_app_pre));
    float *cn_app_cur = (float *)calloc(h_matrix_ref.bits, sizeof(*cn_app_cur));
    float *cn_q_sel_cur = (float *)calloc(h_matrix_ref.bits, sizeof(*cn_q_sel_cur));
    float *cn_r_old_cur = (float *)calloc(h_matrix_ref.bits, sizeof(*cn_r_old_cur));
    float *cn_q_updt_cur = (float *)calloc(h_matrix_ref.bits, sizeof(*cn_q_updt_cur));
    (void)cn_q_sel_cur;

    int **cn_q_sign = (int **)calloc(5 * h_matrix_ref.cols, sizeof(*cn_q_sign));
    for (int i = 0; i < 5 * h_matrix_ref.cols; i++)
        cn_q_sign[i] = (int *)calloc(h_matrix_ref.bits, sizeof(*cn_q_sign[i]));

    int **e_pre = (int **)calloc(h_matrix_ref.rows, sizeof(*e_pre));
    for (int i = 0; i < h_matrix_ref.rows; i++)
        e_pre[i] = (int *)calloc(h_matrix_ref.cols, sizeof(*e_pre[i]));

    ldpc_decoder_output decoder_output;
    printf("use CPU decoder by cdeng\n");

    for (int i = 0; i < h_matrix_ref.rows; i++) {
        for (int j = 0; j < h_matrix_ref.cols; j++)
            e_pre[i][j] = -1;
    }

    for (int i = 0; i < h_matrix_ref.cols; i++) {
        int tmp_pre = -1;
        int first_layer = -1;
        for (int j = 0; j < h_matrix_ref.rows; j++) {
            if ((h_matrix_ref.occupied[j][i] == 1) || ((h_matrix_ref.fade[j][i] == 1) && h_matrix_ref.extra_bytes_of_parity > 0)) {
                e_pre[j][i] = tmp_pre;
                tmp_pre = j;
                if (first_layer < 0)
                    first_layer = j;
            }
        }
        e_pre[first_layer][i] = tmp_pre;
    }

    for (int idx = 0; idx < h_matrix_ref.cols * 8; idx += 1) {
        const auto col_idx = idx / 8;
        const auto word_idx = idx % 8;
        const auto hard_bits = decoder_input.corrupted_codeword.hard.cols[col_idx][word_idx];
        const auto soft_bits_0 = decoder_input.corrupted_codeword.soft0.cols[col_idx][word_idx];
        const auto soft_bits_1 = decoder_input.corrupted_codeword.soft1.cols[col_idx][word_idx];

        for (int bit_idx = 0; bit_idx < 64; bit_idx++) {
            const bool hard = (hard_bits & (1UL << bit_idx));
            const bool sb0 = (soft_bits_0 & (1UL << bit_idx));
            const bool sb1 = (soft_bits_1 & (1UL << bit_idx));
            const auto llr_index = (hard << 2) | (sb1 << 1) | sb0;
            cn_q_mem[col_idx][word_idx * 64 + bit_idx] = ldpc_decoder_parameters.llr_table[llr_index];
            dec_di_blk[col_idx * h_matrix_ref.bits + word_idx * 64 + bit_idx] = llr_index;
        }
    }

    vec_copy(dec_di_blk, dec_do_blk, 0, 0, h_matrix_ref.rows);

    for (int itr = 0; itr < (int)decoder_input.iteration_limit && ((ldec_early_term_en == 0) || (cw_fail == 1)); itr++) {
        cir_cnt = 0;

        for (int layer = 0; layer < h_matrix_ref.rows && ((ldec_early_term_en == 0) || (cw_fail == 1)); layer++) {
            hd_init = (vec_sum(dec_init, h_matrix_ref.cols) != 0);
            cn_c_sel_cur = cn_c_mem[layer];

            for (int i = 0; i < h_matrix_ref.bits; i++) {
                cn_c_updt_cur[i].min1_val = 100000;
                cn_c_updt_cur[i].min2_val = 100000;
                cn_c_updt_cur[i].min1_pos = 0;
                cn_c_updt_cur[i].sign_tot = 1;
            }

            hd_updated = 0;
            vec_clr(layer_synd, h_matrix_ref.bits);

            for (col = 0; col < h_matrix_ref.cols; col++) {
                if (((h_matrix_ref.occupied[layer][col] == 0) && (h_matrix_ref.fade[layer][col] == 0)) || ((h_matrix_ref.fade[layer][col] == 1) && (h_matrix_ref.extra_bytes_of_parity == 0)))
                    continue;

                cn_q_sel_pre = cn_q_mem[col];
                layer_pre = e_pre[layer][col];
                cn_c_sel_pre = cn_c_mem[layer_pre];

                update_node(cn_q_sel_pre, cn_c_sel_pre, cn_r_new_pre, cn_app_pre, layer_pre, col);
                update_hd(dec_init, col, layer, layer_pre, cn_app_pre, cn_app_cur, vn_dec_hd, hd_updated, dec_do_blk, cn_dec_hd, layer_synd);
                update_node_next(cn_r_old_cur, cn_c_sel_cur, cn_q_updt_cur, cn_app_pre, cn_app_cur, ldpc_decoder_parameters.finite_mode, cir_cnt, layer, col, cn_q_sign, cn_q_mem, cn_c_updt_cur);
                cir_cnt++;
            }

            for (int i = 0; i < h_matrix_ref.bits; i++) {
                cn_c_mem[layer][i].min1_val = cn_c_updt_cur[i].min1_val * ldpc_decoder_parameters.alpha;
                cn_c_mem[layer][i].min2_val = cn_c_updt_cur[i].min2_val * ldpc_decoder_parameters.alpha;

                if (ldpc_decoder_parameters.finite_mode == 1) {
                    cn_c_mem[layer][i].min1_val = (float)Sat_Quan((double)cn_c_mem[layer][i].min1_val, ldpc_decoder_parameters.finite_c_max, ldpc_decoder_parameters.finite_c_min, ldpc_decoder_parameters.finite_c_num, ldpc_decoder_parameters.finite_f_num);
                    cn_c_mem[layer][i].min2_val = (float)Sat_Quan((double)cn_c_mem[layer][i].min2_val, ldpc_decoder_parameters.finite_c_max, ldpc_decoder_parameters.finite_c_min, ldpc_decoder_parameters.finite_c_num, ldpc_decoder_parameters.finite_f_num);
                }

                cn_c_mem[layer][i].min1_pos = cn_c_updt_cur[i].min1_pos;
                cn_c_mem[layer][i].sign_tot = cn_c_updt_cur[i].sign_tot;
            }

            layer_synd_wt = vec_sum(layer_synd, h_matrix_ref.bits);
            if (hd_init == 1) {
                hd_stable_cnt = 0;
                synd_pass_cnt = 0;
            } else if ((hd_updated == 0) && (layer_synd_wt == 0)) {
                hd_stable_cnt++;
                synd_pass_cnt++;
            } else {
                hd_stable_cnt = 0;
                synd_pass_cnt = 0;
            }

            if ((synd_pass_cnt >= h_matrix_ref.rows) && (hd_stable_cnt >= (h_matrix_ref.rows - 1))) {
                cw_fail = 0;
                decoder_output.iterations = itr;
            }
        }
    }

    if ((cw_fail == 1) || (ldec_early_term_en == 0))
        decoder_output.iterations = decoder_input.iteration_limit;

    decoder_output.clock_cycles = clock_cycles;
    decoder_output.syndrome_weight_after = layer_synd_wt;
    decoder_output.failure = cw_fail;

    free(dec_init);
    for (int i = 0; i < h_matrix_ref.rows; i++)
        free(cn_c_mem[i]);
    free(cn_c_mem);
    free(cn_c_updt_cur);
    for (int i = 0; i < h_matrix_ref.cols; i++)
        free(cn_q_mem[i]);
    free(cn_q_mem);
    free(cn_r_new_pre);
    free(cn_app_pre);
    free(cn_app_cur);
    free(cn_q_sel_cur);
    free(cn_r_old_cur);
    free(cn_q_updt_cur);
    for (int i = 0; i < 5 * h_matrix_ref.cols; i++)
        free(cn_q_sign[i]);
    free(cn_q_sign);
    for (int i = 0; i < h_matrix_ref.rows; i++)
        free(e_pre[i]);
    free(e_pre);
    free(layer_synd);
    free(cn_dec_hd);
    free(vn_dec_hd);

    printf("decode fail=%d, iter=%d\n", cw_fail, decoder_output.iterations);
    return decoder_output;
}

void logger::log_elapsed_time(const std::chrono::steady_clock::time_point &start_time, const std::chrono::steady_clock::time_point &end_time) {
    static int call_count = 0;
    static std::vector<uint64_t> time_keeper;

    call_count++;
    auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    time_keeper.push_back(elapsed_us);

    if (time_keeper.size() == AVG_WINDOW_SIZE) {
        double sum = 0;
        for (const auto &elem : time_keeper)
            sum += elem;
        std::cout << "  > Decoder Elapsed time: " << (double)sum / time_keeper.size() << " us\n";
        time_keeper.clear();
    }

    return;
}

void logger::print_accumulated_stats(uint64_t accumulated_cw_count, ldpc_decoder_output &decode_stats) {
    double cw_count = accumulated_cw_count;
    std::cout << "===== Codeword Decoding Stats =====\n"
              << " - codewords: " << accumulated_cw_count << "\n"
              << " - failures: " << decode_stats.failure << "\n"
              << " - iterations: " << decode_stats.iterations << "\n"
              << " -> Avg. clock_cycles: " << decode_stats.clock_cycles / cw_count << "\n"
              << " -> Avg. syndrome_weight_before: " << decode_stats.syndrome_weight_before / cw_count << "\n"
              << " -> Avg. syndrome_weight_after: " << decode_stats.syndrome_weight_after / cw_count << "\n"
              << " -> Avg. iterations: " << decode_stats.iterations / cw_count << "\n"
              << " -> Failure rate: " << (double)decode_stats.failure / cw_count << "\n"
              << "\n";
}

void logger::print_accumulated_stats(uint64_t accumulated_cw_count, decoder_output_acc &decode_stats) {
    double cw_count = accumulated_cw_count;
    double total_iter_pass_only = (decode_stats.iterations - logger::ITER_LIMIT * decode_stats.failure);
    double cw_pass_count = (cw_count - decode_stats.failure);

    std::cout << "===== Codeword Decoding Stats =====\n"
              << " - codewords: " << accumulated_cw_count << "\n"
              << " - failures: " << decode_stats.failure << "\n"
              << " - iterations: " << decode_stats.iterations << "\n"
              << " - syndrome_weight_after: " << decode_stats.syndrome_weight_after << "\n"
              << " -> Avg. clock_cycles: " << decode_stats.clock_cycles / cw_count << "\n"
              << " -> Avg. syndrome_weight_before: " << decode_stats.syndrome_weight_before / cw_count << "\n"
              << " -> Avg. syndrome_weight_after: " << decode_stats.syndrome_weight_after / cw_count << "\n"
              << " -> Avg. syndrome_weight_after(Fail only): " << (double)decode_stats.syndrome_weight_after / decode_stats.failure << "\n"
              << " -> Avg. iterations: " << decode_stats.iterations / cw_count << "\n"
              << " -> Avg. iterations (Pass only): " << total_iter_pass_only / cw_pass_count << "\n"
              << " -> Failure rate: " << (double)decode_stats.failure / cw_count << "\n"
              << "\n";
}

void logger::print_iter_stats(uint32_t *iter_info) {
    std::cout << "===== ITER INFO =====\n";
    for (int i = 0, j = 0; i < (int)logger::MAX_ITER; i++) {
        if (iter_info[i] != 0) {
            std::stringstream ss;
            ss << i << ": " << iter_info[i] << ", ";
            std::cout << std::left << std::setw(20) << ss.str();
            j++;
            if (j % 8 == 0)
                std::cout << "\n";
        }
    }
    return;
}
