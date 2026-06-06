#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include "decoder.h"
#include "codec_support.h"
#include "codeword.h"
#include "h_matrix.h"
uint64_t logger::ITER_LIMIT = logger::MAX_ITER;
uint16_t decoder::prng_init[] = {0x9365, 0x49f9, 0xda8f, 0xf7b2, 0x30ee, 0xef08, 0x1b73, 0x8c9a, 0xc646, 0xb550, 0x2edb, 0x71cc, 0x5d27, 0xa8a1, 0x6214, 0x043d, 0x9375, 0x89f9, 0xd98f, 0xf2b2, 0x31ee, 0xef18, 0x0b73, 0x4c9a, 0xc946, 0xba50, 0x3edb, 0x71bc, 0x5327, 0xa8a1, 0x3214, 0x083d};

ldpc_decoder_input::ldpc_decoder_input() {
    return;
}

ldpc_decoder_input::ldpc_decoder_input(uint64_t nand_strobes, uint64_t soft_bits, uint64_t iteration_limit, uint64_t post_iteration) {
    this->nand_strobes = nand_strobes;
    this->soft_bits = soft_bits;
    this->iteration_limit = iteration_limit;
    this->post_iteration = post_iteration;

    return;
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

    return;
}

void decoder_output_acc::print_stats(void) {
    std::cout << "===== Accumulated Decoding Stats =====\n";
    std::cout << "failure: " << failure << "\n";
    std::cout << "iterations: " << iterations << "\n";
    std::cout << "clock_cycles: " << clock_cycles << "\n";
    std::cout << "syndrome_weight_before: " << syndrome_weight_before << "\n";
    std::cout << "syndrome_weight_after: " << syndrome_weight_after << "\n";
    std::cout << "total_errors: " << total_errors << "\n";
    std::cout << "net_acc: " << net_acc << "\n";

    return;
}

decoder::decoder(h_matrix &h_matrix_ref, int nand_strobes, int post_ratio) : h_matrix_ref(h_matrix_ref) {
    ldpc_decoder_parameters.post_ratio = post_ratio;
    ldpc_decoder_parameters.likelihood_init_fraction[0] = 0;
    ldpc_decoder_parameters.likelihood_init_fraction[1] = 0;
    ldpc_decoder_parameters.likelihood_init_fraction[2] = 0;

    if (nand_strobes == 7) {
        ldpc_decoder_parameters.likelihood_init_fraction[1] = 9;
        ldpc_decoder_parameters.likelihood_init_fraction[2] = 3;
    } else {
        ldpc_decoder_parameters.likelihood_init_fraction[1] = 6;
        ldpc_decoder_parameters.likelihood_init_fraction[2] = 6;
    }
}

ldpc_decoder_output decoder::decode_planar(const ldpc_decoder_input &decoder_input) {
    uint64_t iteration = 0;

    uint64_t i, j, k;
    int clock_cycles = 0;
    uint64_t syndrome_weight;
    bool finished = 0;
    bool hamming_weight_lt_circ_thr, hamming_weight_lt_post_thr;
    bool post_trigger1, post_trigger2;

    codeword hard_codeword = decoder_input.corrupted_codeword.hard;

    variable_nodes vn;
    check_nodes cn;
    s_likelihood_levels likelihood_levels;
    ldpc_decoder_output decoder_output;
    s_512_bits prng_512;
    uint64_t syndrome_weight_delayed;
    uint64_t syndrome_weight_r[4] = {0, 0, 0, 0};

    cn = config_check_nodes(h_matrix_ref, hard_codeword);
    syndrome_weight = check_node_weight(cn);

#ifdef ASPEN
    clock_cycles = (h_matrix_ref.cols) + 1;
#else
    clock_cycles = (2 * h_matrix_ref.cols) + 1;
#endif

    if (syndrome_weight == 0)
        finished = 1;

    decoder_output.syndrome_weight_before = syndrome_weight;
    syndrome_weight_delayed = syndrome_weight;

    likelihood_levels = compute_likelihood_levels(decoder_input.nand_strobes, ldpc_decoder_parameters, syndrome_weight, h_matrix_ref.rows);

    for (j = 0; j < h_matrix_ref.cols; j++) {
        for (k = 0; k < h_matrix_ref.bits; k++) {
            int tt = k / 64, ttt = k % 64;
            uint64_t aaa = (decoder_input.corrupted_codeword.soft0.cols[j][tt] >> ttt) & 0x1;
            uint64_t bbb = (decoder_input.corrupted_codeword.soft1.cols[j][tt] >> ttt) & 0x1;
            uint64_t soft_data = aaa | (bbb << 1);

            vn.likelihood[j][k] = likelihood_levels.level[soft_data];
        }
    }

    while (iteration < decoder_input.iteration_limit && !finished) {
        for (j = 0; j < h_matrix_ref.cols; j++) {
            printf("iteration= %llu, col= %llu\n", static_cast<unsigned long long>(iteration), static_cast<unsigned long long>(j));
            clock_cycles++;

            if ((iteration == decoder_input.post_iteration) && (j == 0)) {
                for (i = 0; i < 512; i++) {
                    prng_512.b[i] = (prng_init[int(i / 16)] >> (i % 16)) & 1;
                }
            } else if (iteration >= decoder_input.post_iteration) {
                prng_512 = lfsr_512_bit(prng_512);
            }

            syndrome_weight_r[3] = syndrome_weight_r[2];
            syndrome_weight_r[2] = syndrome_weight_r[1];
            syndrome_weight_r[1] = syndrome_weight_r[0];
            syndrome_weight_r[0] = syndrome_weight;

            syndrome_weight_delayed = (j <= 3) ? syndrome_weight_r[2] : syndrome_weight_r[3];
            syndrome_weight = check_node_weight(cn);

            hamming_weight_lt_circ_thr = (iteration >= decoder_input.post_iteration) && (syndrome_weight_delayed < ldpc_decoder_parameters.syndrome_weight_thr_qc);
            hamming_weight_lt_post_thr = (iteration >= decoder_input.post_iteration) && (syndrome_weight_delayed < ldpc_decoder_parameters.syndrome_weight_thr_post);

            post_trigger1 = hamming_weight_lt_circ_thr && ((iteration % 16) < ldpc_decoder_parameters.post_ratio) && ldpc_decoder_parameters.post_process_en;
            post_trigger2 = hamming_weight_lt_post_thr && ((iteration % 16) >= ldpc_decoder_parameters.post_ratio) && ldpc_decoder_parameters.post_process_en;

            auto start_time = std::chrono::steady_clock::now();

            uint64_t row_mask[LDPC_M][8] = {0};
            uint64_t vn_flip_mask[8];

            for (int word_idx = 0; word_idx < 8; word_idx++) {
                vn_flip_mask[word_idx] = vn.flipped[j][word_idx];
                vn.flipped[j][word_idx] = 0;
            }

            uint64_t active_rows[5] = {LDPC_M, LDPC_M, LDPC_M, LDPC_M, LDPC_M};
            uint64_t active_row_count = 0;

            for (i = 0; i < h_matrix_ref.rows; i++) {
                int idx = h_matrix_ref.operational_h_matrix[i][j];
                h_matrix::range curr_range = h_matrix_ref.tile_ranges[idx];

                if (curr_range.active_count == 0)
                    continue;

                active_rows[active_row_count] = i;
                active_row_count++;

                mask_range_512(row_mask[i], curr_range.offset, curr_range.active_count);
                printf("print row mask[%llu] =", static_cast<unsigned long long>(i));
                for (int ii = 0; ii < 8; ii++)
                    printf(" %016llx ", row_mask[i][7 - ii]);
                printf("\n");
            }

            const uint64_t base_mask_word = 0x0101010101010101ULL;
            const bool base_aggr_condition = (decoder_input.soft_bits > 0) && (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr);

            for (unsigned bit_idx = 0; bit_idx < 8; bit_idx++) {
                const uint64_t active_mask = base_mask_word << bit_idx;
                printf("bit_idx, active_mask=%d   %llx\n", bit_idx, active_mask);

                uint64_t row_mask_m[8];
                uint64_t cn_sample_temp[8];
                uint64_t adder_aligned[8];
                uint64_t weight_word[8] = {0};

                // I. Accumulate weight from all check_node rows
                for (uint64_t row_num = 0; row_num < active_row_count; row_num++) {
                    i = active_rows[row_num];
                    // 1) Limit the mask to active bit locations of this run.
                    for (int word_idx = 0; word_idx < 8; word_idx++)
                        cn_sample_temp[word_idx] = row_mask[i][word_idx] & active_mask;

                    printf("print cn_sample_temp[%llu] =", static_cast<unsigned long long>(i));
                    for (int ii = 0; ii < 8; ii++)
                        printf(" %016llx ", cn_sample_temp[7 - ii]);
                    printf("\n");
                    // 2) shift by m offset
                    rotate_left_512(cn_sample_temp, row_mask_m, h_matrix_ref.element[i][j]);

                    // 3) sample cn[i]
                    for (int word_idx = 0; word_idx < 8; word_idx++)
                        cn_sample_temp[word_idx] = row_mask_m[word_idx] & cn.rows[i][word_idx];

                    // 4) shift back to align weight sums. Have additional left shift so we can align all with adder base.
                    rotate_left_512(cn_sample_temp, adder_aligned, -h_matrix_ref.element[i][j] + bit_idx);

                    // 5) Accumulate results into weights
                    for (int word_idx = 0; word_idx < 8; word_idx++)
                        weight_word[word_idx] += adder_aligned[word_idx];
                }

                if (base_aggr_condition) {
                    for (int word_idx = 0; word_idx < 8; word_idx++) {
                        uint64_t curr_aggr = vn_flip_mask[word_idx] & active_mask;
                        uint8_t *aggr_ptr = (uint8_t *)&curr_aggr;
                        uint8_t *weight_ptr = (uint8_t *)&weight_word[word_idx];

                        for (int byte_idx = 0; byte_idx < 8; byte_idx++) {
                            const bool update_weight = (aggr_ptr[byte_idx] != 0) && (weight_ptr[byte_idx] > 0);
                            weight_ptr[byte_idx] += update_weight * (weight_ptr[byte_idx] - 1);
                        }
                    }
                }

                update_vn_planar(vn.flipped[j], (uint8_t *)weight_word, &vn.likelihood[j][0], likelihood_levels, prng_512.b, post_trigger1, post_trigger2, bit_idx);
            }

            for (int word_idx = 0; word_idx < 8; word_idx++)
                vn_flip_mask[word_idx] ^= vn.flipped[j][word_idx];

            for (uint64_t row_num = 0; row_num < active_row_count; row_num++) {
                i = active_rows[row_num];
                uint64_t vn_flip_mask_m_aligned[8];

                for (int word_idx = 0; word_idx < 8; word_idx++)
                    row_mask[i][word_idx] &= vn_flip_mask[word_idx];

                rotate_left_512(row_mask[i], vn_flip_mask_m_aligned, h_matrix_ref.element[i][j]);

                for (int word_idx = 0; word_idx < 8; word_idx++)
                    cn.rows[i][word_idx] ^= vn_flip_mask_m_aligned[word_idx];
            }

            auto end_time = std::chrono::steady_clock::now();
        }

        finished = (check_node_weight(cn) == 0);
        clock_cycles++;
        iteration++;
    }

    decoder_output.iterations = iteration;
    decoder_output.clock_cycles = clock_cycles;
    decoder_output.syndrome_weight_after = check_node_weight(cn);
    decoder_output.failure = (decoder_output.syndrome_weight_after != 0);

    for (j = 0; j < h_matrix_ref.cols; j++) {
        for (int word_idx = 0; word_idx < 8; word_idx++)
            decoder_output.total_errors += __builtin_popcountll(vn.flipped[j][word_idx]);
    }

    return decoder_output;
}

void decoder::update_vn_planar(uint64_t *vn_flipped, const uint8_t *weight_word, unsigned char *vn_likelihood, const s_likelihood_levels &likelihood_levels, const bool *prng_512, bool post_trigger1, bool post_trigger2, int bit_offset) {
    for (int weight_idx = 0; weight_idx < 64; weight_idx++) {
        const auto k = weight_idx * 8 + bit_offset;
        const auto curr_likelihood = vn_likelihood[k];
        const bool flipped_prev = (curr_likelihood >= likelihood_levels.flip_thr);
        int likelihood_new = flipped_prev ? curr_likelihood - weight_word[weight_idx] : curr_likelihood + weight_word[weight_idx] - 1;

        if (post_trigger1 && prng_512[k] && (likelihood_new <= likelihood_levels.flip_thr)) {
            likelihood_new = likelihood_levels.flip_thr - (likelihood_new < likelihood_levels.flip_thr) + (likelihood_new == likelihood_levels.flip_thr);
        }

        likelihood_new += post_trigger2 && prng_512[k] && !flipped_prev && (weight_word[weight_idx] == 1) && (likelihood_new == (likelihood_levels.flip_thr - 1));

        if (likelihood_new <= likelihood_levels.min)
            likelihood_new = likelihood_levels.min;
        else if (likelihood_new >= likelihood_levels.max)
            likelihood_new = likelihood_levels.max;

        vn_likelihood[k] = likelihood_new;

        const uint64_t new_active = (likelihood_new >= likelihood_levels.flip_thr) ? UNIT : 0ULL;
        vn_flipped[k / 64] |= new_active << (k % 64);
    }
}

s_likelihood_levels decoder::compute_likelihood_levels(int strobes, decoder::ldpc_decoder_params ldpc_decoder_parameters, int syndrome_weight, int rows) {
    s_likelihood_levels likelihood_levels;
    int address;
    int delta[4];
    int delta_sum;
    int delta_total;
    int coef[4];
    int weak_minus_strong;
    const unsigned VN_BITS = decoder::ldpc_decoder_params::VN_BITS;

    likelihood_levels.max = (1U << VN_BITS) - 1;
    likelihood_levels.min = 1;
    likelihood_levels.flip_thr = (VN_BITS == 3) ? likelihood_levels.max - 3 : likelihood_levels.max - 7;
    likelihood_levels.weak = likelihood_levels.flip_thr - 4;
    likelihood_levels.strong = likelihood_levels.weak;
    likelihood_levels.level[0] = likelihood_levels.weak;
    likelihood_levels.level[1] = likelihood_levels.weak;
    likelihood_levels.level[2] = likelihood_levels.weak;
    likelihood_levels.level[3] = likelihood_levels.weak;

    if (strobes > 0) {
        address = syndrome_weight >> 5;
        if (address >= 64)
            address = 63;

        int coef_index = 0;
        if (rows == 7)
            coef_index = 1;
        if (rows == 8)
            coef_index = 2;
        if (rows == 9)
            coef_index = 3;
        if (rows == 10)
            coef_index = 4;
        if (rows == 11)
            coef_index = 5;
        if (rows == 12)
            coef_index = 6;
        if (rows == 13)
            coef_index = 7;

        coef[0] = ldpc_decoder_parameters.likelihood_init_coef_all[coef_index][0];
        coef[1] = ldpc_decoder_parameters.likelihood_init_coef_all[coef_index][1];
        coef[2] = ldpc_decoder_parameters.likelihood_init_coef_all[coef_index][2];
        coef[3] = ldpc_decoder_parameters.likelihood_init_coef_all[coef_index][3];

        delta[0] = (address >= 0) ? (coef[0] * (address - 0)) : 0;
        delta[1] = (address >= 8) ? (coef[1] * (address - 8)) : 0;
        delta[2] = (address >= 16) ? (coef[2] * (address - 16)) : 0;
        delta[3] = (address >= 24) ? (coef[3] * (address - 24)) : 0;

        delta_sum = delta[0] + delta[1] + delta[2] + delta[3];
        delta_total = (VN_BITS == 8) ? ((delta_sum >> 1) + (delta_sum >> 2)) : ((delta_sum >> 2) + (delta_sum >> 3));
        likelihood_levels.strong = likelihood_levels.weak - delta_total;

        if (likelihood_levels.strong < likelihood_levels.min)
            likelihood_levels.strong = likelihood_levels.min;

        likelihood_levels.level[0] = likelihood_levels.strong;
        likelihood_levels.level[1] = likelihood_levels.level[0];
        likelihood_levels.level[2] = likelihood_levels.level[3];

        if (strobes > 3) {
            likelihood_levels.level[1] = likelihood_levels.level[0] + ldpc_decoder_parameters.likelihood_init_fraction[1];
            likelihood_levels.level[2] = likelihood_levels.level[3] - ldpc_decoder_parameters.likelihood_init_fraction[2];
            weak_minus_strong = likelihood_levels.level[3] - likelihood_levels.level[0];

            if (weak_minus_strong >= 8)
                likelihood_levels.min = likelihood_levels.strong;
        }
    }

    std::cout << ":: likelihood levels"
              << "\n"
              << "level[0]:" << likelihood_levels.level[0] << "\n"
              << "level[1]:" << likelihood_levels.level[1] << "\n"
              << "level[2]:" << likelihood_levels.level[2] << "\n"
              << "level[3]:" << likelihood_levels.level[3] << "\n"
              << "Flip thr: " << likelihood_levels.flip_thr << "\n"
              << "\n";

    return likelihood_levels;
}

s_512_bits decoder::lfsr_512_bit(s_512_bits data_in) {
    s_512_bits data_out;
    int i;
    bool bit0;

    bit0 = data_in.b[0];
    for (i = 0; i < 512; i++) {
        if (i == 511)
            data_out.b[i] = bit0;
        else if ((i == 509) || (i == 506) || (i == 503))
            data_out.b[i] = data_in.b[i + 1] ^ bit0;
        else
            data_out.b[i] = data_in.b[i + 1];
    }

    return data_out;
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
        std::cout << " -> Decoder Elapsed time: " << (double)sum / time_keeper.size() << " us\n";
        time_keeper.clear();
    }

    return;
}

void logger::print_accumulated_stats(uint64_t accumulated_cw_count, ldpc_decoder_output &decode_stats) {
    double cw_count = accumulated_cw_count;

    if (cw_count == 0)
        cw_count = 1;

    std::cout << "===== Codeword Decoding Stats =====\n"
              << " - codewords:" << accumulated_cw_count << "\n"
              << " - failures:" << decode_stats.failure << "\n"
              << " - iterations:" << decode_stats.iterations << "\n"
              << " - Avg.clock_cycles:" << decode_stats.clock_cycles / cw_count << "\n"
              << " ->Avg. syndrome_weight_before: " << decode_stats.syndrome_weight_before / cw_count << "\n"
              << " ->Avg. syndrome_weight_after:" << decode_stats.syndrome_weight_after / cw_count << "\n"
              << " -> Avg. iterations: " << decode_stats.iterations / cw_count << "\n"
              << " -> Failure rate: " << (double)decode_stats.failure / cw_count << "\n"
              << "\n";
}

void logger::print_accumulated_stats(uint64_t accumulated_cw_count, decoder_output_acc &decode_stats) {
    double cw_count = accumulated_cw_count;

    if (cw_count == 0)
        cw_count = 1;

    double cw_pass_count = cw_count - decode_stats.failure;
    double total_iter_pass_only = decode_stats.iterations - logger::ITER_LIMIT * decode_stats.failure;
    double fail_count = decode_stats.failure;

    if (cw_pass_count == 0)
        cw_pass_count = 1;
    if (fail_count == 0)
        fail_count = 1;

    std::cout << "==== Codeword Decoding Stats ====\n"
              << " - codewords:" << accumulated_cw_count << "\n"
              << " - failures:" << decode_stats.failure << "\n"
              << " - iterations:" << decode_stats.iterations << "\n"
              << " - sum syndrome_weight_after:" << decode_stats.syndrome_weight_after << "\n"
              << " - Avg.clock_cycles:" << decode_stats.clock_cycles / cw_count << "\n"
              << " ->Avg. syndrome_weight_before: " << decode_stats.syndrome_weight_before / cw_count << "\n"
              << " ->Avg. syndrome_weight_after: " << decode_stats.syndrome_weight_after / cw_count << "\n"
              << " -> Avg. syndrome_weight_after (Fail only): " << (double)decode_stats.syndrome_weight_after / fail_count << "\n"
              << " -> Avg. iterations:" << decode_stats.iterations / cw_count << "\n"
              << " -> Avg. iterations (Pass only): " << total_iter_pass_only / cw_pass_count << "\n"
              << " -> Failure rate: " << (double)decode_stats.failure / cw_count << "\n"
              << "\n";
}

void logger::print_iter_stats(uint32_t *iter_info) {
    std::cout << "\n-------------------------------- ITER INFO ------------------------------\n";

    for (int i = 0, j = 0; i < logger::MAX_ITER; i++) {
        if (iter_info[i] != 0) {
            std::stringstream ss;
            ss << i << ":" << iter_info[i] << ",";
            std::cout << std::left << std::setw(20) << ss.str();
            j++;

            if (j % 8 == 0)
                std::cout << "\n";
        }
    }

    return;
}
