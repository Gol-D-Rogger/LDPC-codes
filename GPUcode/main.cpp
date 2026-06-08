#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "codeword.h"
#include "h_matrix.h"
#include "error_injector.h"
#include "decoder.h"
#include "decoder_gpu.h"

void ch_llr_alloc(int sd_bit, const float *vref_in, int llr_tot_bit, int llr_frac_bit, float hd0_llr, float hd1_llr, float rber, float llr_table[8], float &awgn_sigma);

struct program_config {
    uint64_t total_codewords = 5000000;
    uint64_t max_failure_count = 100;
    int bytes_of_userdata = 4112;
    int bytes_of_parity = 472;
    int nand_strobes = 7;
    int err_inj_mode = 0;
    int iteration_limit = 32;
    int post_iteration = 10;
    int post_ratio = 0;
    int mx_cnfg = 0;
    int finite_mode = 1;
    int finite_q_num = 9;
    int finite_f_num = 3;
    int finite_c_num = 7;
    int llr_tot_bit = 7;
    int llr_frac_bit = 3;
    float vref[7] = {0, 0.15, -0.15, 0.3, -0.3, 0.5, -0.5};
    float hd0_llr = 1.875; float hd1_llr = -1.875;
    double rber = 0.015;
    double target_distribution_rber = 0.025;
    bool save_decode_fail_data = false;
    bool GPU_mode = true;
    float alpha = 0.625;

    int parse_args(int argc, char **argv);
};

program_config Config;

int program_config::parse_args(int argc, char **argv) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--", 0) == 0) {
            std::string key = arg.substr(2);

            if (i + 1 < argc) {
                std::string value = argv[++i];

                if (key == "total_codewords")
                    total_codewords = std::stoull(value);
                else if (key == "max_failure_count")
                    max_failure_count = std::stoull(value);
                else if (key == "bytes_of_userdata")
                    bytes_of_userdata = std::stoi(value);
                else if (key == "bytes_of_parity")
                    bytes_of_parity = std::stoi(value);
                else if (key == "nand_strobes")
                    nand_strobes = std::stoi(value);
                else if (key == "err_inj_mode")
                    err_inj_mode = std::stoi(value);
                else if (key == "iteration_limit")
                    iteration_limit = std::stoi(value);
                // else if (key == "post_iteration")
                //     post_iteration = std::stoi(value);
                // else if (key == "post_ratio")
                //     post_ratio = std::stoi(value);
                else if (key == "mx_cnfg")
                    mx_cnfg = std::stoi(value);
                else if (key == "rber")
                    rber = std::stod(value);
                else if (key == "target_distribution_rber")
                    target_distribution_rber = std::stod(value);
                else if (key == "save_decode_fail_data")
                    save_decode_fail_data = std::stoi(value);
                else if (key == "GPU_mode")
                    GPU_mode = std::stoi(value);
                else if (key == "finite_mode")
                    finite_mode = std::stoi(value);
                else if (key == "finite_q_num")
                    finite_q_num = std::stoi(value);
                else if (key == "finite_f_num")
                    finite_f_num = std::stoi(value);
                else if (key == "finite_c_num")
                    finite_c_num = std::stoi(value);
                else if (key == "llr_tot_bit")
                    llr_tot_bit = std::stoi(value);
                else if (key == "llr_frac_bit")
                    llr_frac_bit = std::stoi(value);
                else {
                    std::cerr << "Unknown argument: " << key << std::endl;
                    exit(0);
                }
            }
        }
    }

    return 0;
}

class configurator {
  public:
    const program_config &config_ref;
    uint64_t batch_codewords = 1;
    int base_matrix_rows;
    int base_matrix_cols;
    int soft_bits;
    h_matrix *h_matrix_ptr;
    error_injector *error_injector_ptr;
    decoder *decoder_ptr;
    ldpc_decoder_input *ldpc_decoder_input_ptr;

    configurator(const program_config &config = Config) : config_ref(config) {
        if (config.nand_strobes >= 5)
            soft_bits = 2;
        else if (config.nand_strobes >= 3)
            soft_bits = 1;
        else
            soft_bits = 0;

        base_matrix_rows = config.bytes_of_parity / 64 + ((config.bytes_of_parity % 64) != 0);
        base_matrix_cols = config.bytes_of_userdata / 64 + ((config.bytes_of_userdata % 64) != 0);
        base_matrix_cols = base_matrix_rows + base_matrix_cols;

        h_matrix_ptr = new h_matrix(config.bytes_of_userdata, config.bytes_of_parity, config.mx_cnfg);

        int parity_col_idx = base_matrix_cols - base_matrix_rows;
        error_injector_ptr = new error_injector(config.err_inj_mode, config.nand_strobes, config.rber, soft_bits, parity_col_idx, h_matrix_ptr->extra_bits_of_userdata, h_matrix_ptr->extra_bits_of_parity, base_matrix_cols, config.target_distribution_rber);

        std::printf("initialize decoder\n");
        decoder_ptr = new decoder(*h_matrix_ptr, config.nand_strobes, config.post_ratio);

        ldpc_decoder_input_ptr = new ldpc_decoder_input(config.nand_strobes, soft_bits, config.iteration_limit, config.post_iteration);

        if (config_ref.GPU_mode)
            configure_gpu_info();

        std::cout << "[CONFIGURATOR] :: Configuration complete." << std::endl;
        return;
    }

    ~configurator() {
        std::cout << "[CONFIGURATOR] :: Releasing resources." << std::endl;
        delete decoder_ptr;
        delete h_matrix_ptr;
        delete error_injector_ptr;
        delete ldpc_decoder_input_ptr;
    }

    void print_configurations(void) {
        std::cout << "\n==== Program Parameters Configured ====" << "\n";
        std::cout << " - total_codewords: " << config_ref.total_codewords << "\n";
        std::cout << " - bytes_of_userdata: " << config_ref.bytes_of_userdata << "\n";
        std::cout << " - bytes_of_parity: " << config_ref.bytes_of_parity << "\n";
        std::cout << " - nand_strobes: " << config_ref.nand_strobes << "\n";
        std::cout << " - max_failure_count: " << config_ref.max_failure_count << "\n";
        std::cout << " - rber: " << config_ref.rber << "\n";
        std::cout << " - save_decode_fail_data: " << config_ref.save_decode_fail_data << "\n";
        std::cout << " -> base_matrix_rows: " << base_matrix_rows << "\n";
        std::cout << " -> base_matrix_cols: " << base_matrix_cols << "\n";
        std::cout << " -> soft_bits: " << soft_bits << "\n";
        std::cout << " - GPU_mode: " << config_ref.GPU_mode << "\n";
        std::cout << std::endl;
    }

  private:
    void configure_gpu_info(void) {
        device_global_info_struct host_decoder_info_struct{};

        host_decoder_info_struct.cols = h_matrix_ptr->cols;
        host_decoder_info_struct.rows = h_matrix_ptr->rows;

        host_decoder_info_struct.nand_strobes = ldpc_decoder_input_ptr->nand_strobes;
        host_decoder_info_struct.soft_bits = ldpc_decoder_input_ptr->soft_bits;
        host_decoder_info_struct.post_iteration = ldpc_decoder_input_ptr->post_iteration;
        host_decoder_info_struct.iteration_limit = ldpc_decoder_input_ptr->iteration_limit;
        std::printf("host_decoder_info_struct.iteration_limit = %d\n", host_decoder_info_struct.iteration_limit);

        host_decoder_info_struct.VN_BITS = decoder_ptr->ldpc_decoder_parameters.VN_BITS;
        host_decoder_info_struct.post_process_en = decoder_ptr->ldpc_decoder_parameters.post_process_en;
        // host_decoder_info_struct.syndrome_weight_thr_qc = decoder_ptr->ldpc_decoder_parameters.syndrome_weight_thr_qc;
        // host_decoder_info_struct.syndrome_weight_thr_post = decoder_ptr->ldpc_decoder_parameters.syndrome_weight_thr_post;
        // host_decoder_info_struct.likelihood_thr = decoder_ptr->ldpc_decoder_parameters.likelihood_thr;
        // host_decoder_info_struct.post_ratio = decoder_ptr->ldpc_decoder_parameters.post_ratio;
        // std::memcpy(host_decoder_info_struct.likelihood_init_coef_all, decoder_ptr->ldpc_decoder_parameters.likelihood_init_coef_all, sizeof(host_decoder_info_struct.likelihood_init_coef_all));
        // std::memcpy(host_decoder_info_struct.likelihood_init_fraction, decoder_ptr->ldpc_decoder_parameters.likelihood_init_fraction, sizeof(host_decoder_info_struct.likelihood_init_fraction));

        for (int i = 0; i < error_injector::MAX_THRESHOLDS; i++)
            host_decoder_info_struct.error_region_prob[i] = error_injector_ptr->error_region_prob[i];

        host_decoder_info_struct.max = (1ULL << host_decoder_info_struct.VN_BITS) - 1;
        host_decoder_info_struct.flip_thr = (host_decoder_info_struct.VN_BITS == 3) ? host_decoder_info_struct.max - 3 : host_decoder_info_struct.max - 7;

        for (int i = 0; i < h_matrix_ptr->rows; i++) {
            for (int j = 0; j < h_matrix_ptr->cols; j++)
                host_decoder_info_struct.e_pre[i][j] = -1;
        }

        for (int i = 0; i < h_matrix_ptr->cols; i++) {
            int tmp_pre = -1;
            int first_layer = -1;
            for (int j = 0; j < h_matrix_ptr->rows; j++) {
                if ((h_matrix_ptr->occupied[j][i] == 1) || ((h_matrix_ptr->fade[j][i] == 1) && h_matrix_ptr->extra_bytes_of_parity > 0)) {
                    host_decoder_info_struct.e_pre[j][i] = tmp_pre;
                    tmp_pre = j;
                    if (first_layer < 0)
                        first_layer = j;
                }
            }

            if (first_layer >= 0)
                host_decoder_info_struct.e_pre[first_layer][i] = tmp_pre;
        }

        host_decoder_info_struct.bits = h_matrix_ptr->bits;
        host_decoder_info_struct.extra_bytes_of_parity = h_matrix_ptr->extra_bytes_of_parity;
        host_decoder_info_struct.extra_bits_of_parity = h_matrix_ptr->extra_bits_of_parity;
        host_decoder_info_struct.extra_bytes_of_userdata = h_matrix_ptr->extra_bytes_of_userdata;

        for (int i = 0; i < h_matrix_ptr->rows; i++) {
            for (int j = 0; j < h_matrix_ptr->cols; j++) {
                host_decoder_info_struct.element[i][j] = (uint16_t)h_matrix_ptr->element[i][j];
                host_decoder_info_struct.fade[i][j] = (char)h_matrix_ptr->fade[i][j];
                host_decoder_info_struct.occupied[i][j] = (char)h_matrix_ptr->occupied[i][j];
            }
        }

        for (int col = 0; col < h_matrix_ptr->cols; col++) {
            for (int i = 0; i < 32; ++i) {
                uint16_t packed = 0;
                for (int bit = 0; bit < 16; ++bit)
                    packed |= (static_cast<uint16_t>(h_matrix_ptr->last_row_active_bits[col][i * 16 + bit]) << bit);
                host_decoder_info_struct.mask[col][i] = packed;
            }
        }

        for (int i = 0; i < 7; i++)
            host_decoder_info_struct.vref[i] = config_ref.vref[i];

        host_decoder_info_struct.alpha = config_ref.alpha;
        host_decoder_info_struct.finite_mode = config_ref.finite_mode;
        host_decoder_info_struct.finite_q_num = config_ref.finite_q_num;
        host_decoder_info_struct.finite_c_num = config_ref.finite_c_num;
        host_decoder_info_struct.finite_f_num = config_ref.finite_f_num;

        float tmp = std::pow(2.0f, config_ref.finite_f_num);
        host_decoder_info_struct.finite_q_max = (std::pow(2.0f, config_ref.finite_q_num - 1) - 1) / tmp;
        host_decoder_info_struct.finite_q_min = (1 - std::pow(2.0f, config_ref.finite_q_num - 1)) / tmp;
        host_decoder_info_struct.finite_r_max = (std::pow(2.0f, config_ref.finite_c_num - 1) - 1) / tmp;
        host_decoder_info_struct.finite_r_min = (1 - std::pow(2.0f, config_ref.finite_c_num - 1)) / tmp;
        host_decoder_info_struct.finite_c_max = (std::pow(2.0f, config_ref.finite_c_num - 1) - 1) / tmp;
        host_decoder_info_struct.finite_c_min = 0;

        std::printf("finite_q_num = %d, finite_c_num=%d, finite_f_num=%d,llr_tot_bit=%d,llr_frac_bit=%d\n", config_ref.finite_q_num, config_ref.finite_c_num, config_ref.finite_f_num, config_ref.llr_tot_bit, config_ref.llr_frac_bit);
        std::printf("q max/min = %f %f\n", host_decoder_info_struct.finite_q_max, host_decoder_info_struct.finite_q_min);
        std::printf("r max/min = %f %f\n", host_decoder_info_struct.finite_r_max, host_decoder_info_struct.finite_r_min);
        std::printf("c max/min = %f %f\n", host_decoder_info_struct.finite_c_max, host_decoder_info_struct.finite_c_min);

        ch_llr_alloc(config_ref.nand_strobes, &config_ref.vref[0], config_ref.llr_tot_bit, config_ref.llr_frac_bit, config_ref.hd0_llr, config_ref.hd1_llr, config_ref.rber, host_decoder_info_struct.llr_table, host_decoder_info_struct.awgn_sigma);
        std::printf("awgn_sigma = %f \n", host_decoder_info_struct.awgn_sigma);

        copy_decoder_info_to_device(host_decoder_info_struct);
        std::cout << "[CONFIGURATOR] : Host side device_global_info_struct configured." << std::endl;
    }
};

uint64_t test_codewords_batch(configurator &configurator_obj, uint64_t &accumulated_cw_count, ldpc_decoder_output &accumulated_decode_stats) {
    const program_config &cfg = configurator_obj.config_ref;
    error_injector &error_injector_ref = *configurator_obj.error_injector_ptr;
    decoder &decoder_ref = *configurator_obj.decoder_ptr;
    ldpc_decoder_input &decoder_input = *configurator_obj.ldpc_decoder_input_ptr;
    uint64_t curr_batch_cw_idx;

    for (curr_batch_cw_idx = 0; curr_batch_cw_idx < configurator_obj.batch_codewords; curr_batch_cw_idx++) {
        decoder_input.clear_cw();
        error_injector_ref.inject_errors(decoder_input.corrupted_codeword);
        ldpc_decoder_output decoder_output = decoder_ref.decode_planar(decoder_input);
        accumulated_decode_stats.failure += decoder_output.failure;
        accumulated_decode_stats.total_errors += decoder_output.total_errors;
        accumulated_decode_stats.iterations += decoder_output.iterations;
        accumulated_decode_stats.clock_cycles += decoder_output.clock_cycles;
        accumulated_decode_stats.syndrome_weight_before += decoder_output.syndrome_weight_before;
        accumulated_decode_stats.syndrome_weight_after += decoder_output.syndrome_weight_after;

        if (accumulated_decode_stats.failure > cfg.max_failure_count) {
            curr_batch_cw_idx += 1;
            break;
        }
    }

    accumulated_cw_count += curr_batch_cw_idx;
    return curr_batch_cw_idx;
}

int cpu_main(configurator &configurator_obj) {
    uint64_t accumulated_cw_count = 0;
    ldpc_decoder_output accumulated_decode_stats;

    while (accumulated_cw_count < Config.total_codewords) {
        auto start_time = std::chrono::steady_clock::now();
        uint64_t current_batch_count = test_codewords_batch(configurator_obj, accumulated_cw_count, accumulated_decode_stats);
        auto end_time = std::chrono::steady_clock::now();

        logger::print_accumulated_stats(accumulated_cw_count, accumulated_decode_stats);

        if (current_batch_count > 0) {
            std::chrono::duration<double, std::milli> duration = end_time - start_time;
            double time_per_codeword_ms = duration.count() / current_batch_count;
            std::cout << " -> Time per codeword(ms): " << std::fixed << std::setprecision(6) << time_per_codeword_ms << "\n";
        }

        if (accumulated_decode_stats.failure > Config.max_failure_count)
            break;
    }

    std::cout << "Tested " << accumulated_cw_count << " codewords\n";
    return 0;
}

int main(int argc, char **argv) {
    Config.parse_args(argc, argv);
    configurator configurator_obj;
    configurator_obj.print_configurations();

    if (LDPC_N > 96) {
        std::printf("Error, LDPC_N too large, need to change dec_init[3] in OptimizedSharedMemory\n");
        return 1;
    }

    if (!Config.GPU_mode)
        return cpu_main(configurator_obj);

    error_injector &error_injector_ref = *configurator_obj.error_injector_ptr;
    decoder &decoder_ref = *configurator_obj.decoder_ptr;
    ldpc_decoder_input &decoder_input = *configurator_obj.ldpc_decoder_input_ptr;

    (void)decoder_ref;
    std::printf("Error injector Seed: %lu\n", error_injector_ref.seed);

    decoder_input.clear_cw();
    error_injector_ref.inject_errors(decoder_input.corrupted_codeword);

    std::printf("call before execute_gpu_kernel\n");
    execute_gpu_kernel(Config.total_codewords, Config.max_failure_count);

    return 0;
}
