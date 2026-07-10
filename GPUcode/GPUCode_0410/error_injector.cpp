#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "error_injector.h"

error_injector::error_injector() {
    seed = RANDOM_SEED ? rd() : seed;
    gen.seed(seed);
    return;
}

error_injector::error_injector(int insertion_mode, int strobes, double rber, int soft_bits, int parity_col_idx, int final_udata_col_bits, int final_parity_col_bits, int codeword_size, int bytes_of_userdata, int bytes_of_parity, double target_distribution_rber) {
    this->insertion_mode = insertion_mode;
    this->strobes = strobes;
    this->rber = rber;
    this->parity_col_idx = parity_col_idx;
    this->final_udata_col_bits = final_udata_col_bits;
    this->final_parity_col_bits = final_parity_col_bits;
    this->codeword_size = codeword_size;
    this->bytes_of_userdata = bytes_of_userdata;
    this->bytes_of_parity = bytes_of_parity;
    this->target_distribution_rber = target_distribution_rber;
    this->num_soft_bits = soft_bits;

    seed = RANDOM_SEED ? rd() : seed;
    gen.seed(seed);

    std::cout << "[ERROR_INJECTOR] :: Seed: " << seed << "\n";

    if (insertion_mode == 0) {
        inject_errors = std::bind(&error_injector::inject_normal_distribution_errors, this, std::placeholders::_1);
        std::cout << "[ERROR_INJECTOR] :: Normal Error Injection mode configured.\n";
    } else {
        inject_errors = std::bind(&error_injector::inject_fixed_count_errors, this, std::placeholders::_1);
        std::cout << "[ERROR_INJECTOR] :: Fixed count Error Injection mode configured.\n";
    }

    if (prob_region_mode == 0)
        set_error_region_probs_tbench(rber);
    else
        set_error_region_probs_cmodel(rber);

    std::cout << "[ERROR_INJECTOR] :: " << (prob_region_mode == 0 ? "TEST BENCH" : "C-MODEL") << " error regioning chosen.\n";
    std::cout << "[ERROR_INJECTOR] :: Error probabilities set with " << strobes << " strobes and RBER = " << rber << ".\n";
    std::cout << "[ERROR_INJECTOR] :: Target distribution RBER = " << target_distribution_rber << "\n";
    std::cout << "[ERROR_INJECTOR] :: Error probabilities:\n";
    for (int i = 0; i < MAX_THRESHOLDS; i++)
        std::cout << "  error_region_prob[" << i << "] = " << error_region_prob[i] << "\n";

    if (strobes == 3) {
        std::cout << "[ERROR_INJECTOR] :: Symbol region probability (Regioned for strobe 3)\n"
                  << "  probability_strong_correct = " << 1 - error_region_prob[1] << "\n"
                  << "  probability_weak_correct = " << error_region_prob[1] - error_region_prob[0] << "\n"
                  << "  probability_weak_error = " << error_region_prob[0] - error_region_prob[2] << "\n"
                  << "  probability_strong_error = " << error_region_prob[2] << "\n"
                  << "  Net_error_probab(RBER) = " << error_region_prob[0] << "\n";
    } else if (strobes == 7) {
        std::cout << "[ERROR_INJECTOR] :: Regioning Thresholds (Regioned for strobe 7)\n"
                  << "  Center = " << thresholds[0] - thresholds[0] << "\n"
                  << "  +1st = " << thresholds[2] - thresholds[0] << "\n"
                  << "  +2nd = " << thresholds[4] - thresholds[0] << "\n"
                  << "  +3rd = " << thresholds[6] - thresholds[0] << "\n"
                  << "  -1st = " << thresholds[1] - thresholds[0] << "\n"
                  << "  -2nd = " << thresholds[3] - thresholds[0] << "\n"
                  << "  -3rd = " << thresholds[5] - thresholds[0] << "\n";
        std::cout << "[ERROR_INJECTOR] :: Symbol region probability (Regioned for strobe 7)\n"
                  << "  P(symbol correct; strength 0) = " << error_region_prob[1] - error_region_prob[0] << "\n"
                  << "  P(symbol correct; strength 1) = " << error_region_prob[3] - error_region_prob[1] << "\n"
                  << "  P(symbol correct; strength 2) = " << error_region_prob[5] - error_region_prob[3] << "\n"
                  << "  P(symbol correct; strength 3) = " << 1 - error_region_prob[5] << "\n"
                  << "  P(symbol wrong ; strength 0) = " << error_region_prob[0] - error_region_prob[2] << "\n"
                  << "  P(symbol wrong ; strength 1) = " << error_region_prob[2] - error_region_prob[4] << "\n"
                  << "  P(symbol wrong ; strength 2) = " << error_region_prob[4] - error_region_prob[6] << "\n"
                  << "  P(symbol wrong ; strength 3) = " << error_region_prob[6] << "\n";
    }

    std::cout << "\n";
    return;
}

bool error_injector::bit_active(int col_idx, int bit_idx) const {
    const int bits_per_col = codeword::WORD_COUNT * codeword::WORD_SIZE;
    if (col_idx < parity_col_idx)
        return (col_idx * bits_per_col + bit_idx) < (bytes_of_userdata * 8);

    return ((col_idx - parity_col_idx) * bits_per_col + bit_idx) < (bytes_of_parity * 8);
}

void error_injector::inject_normal_distribution_errors(decoder_input_cw &cw) {
    const uint64_t UNIT = 1;
    std::uniform_real_distribution<double> dist(0.0f, 1.0f);

    for (int col_idx = 0; col_idx < codeword_size; col_idx++) {
        for (int word_idx = 0; word_idx < codeword::WORD_COUNT; word_idx++) {
            uint64_t error_pattern = 0;
            uint64_t soft0_pattern = 0;
            uint64_t soft1_pattern = 0;

            for (unsigned bit_idx = 0; bit_idx < 64; bit_idx++) {
                int k = word_idx * 64 + bit_idx;
                bool blank_bit = !bit_active(col_idx, k);

                double vn = blank_bit ? 1.00 : dist(gen);

                if (vn < rber)
                    error_pattern |= UNIT << bit_idx;

                if (num_soft_bits == 1) {
                    if ((this->error_region_prob[2] < vn) && (vn < this->error_region_prob[1])) {
                        soft0_pattern |= UNIT << bit_idx;
                        soft1_pattern |= UNIT << bit_idx;
                    }
                } else if (num_soft_bits == 2) {
                    if ((this->error_region_prob[2] < vn) && (vn < this->error_region_prob[1])) {
                        soft0_pattern |= UNIT << bit_idx;
                        soft1_pattern |= UNIT << bit_idx;
                    } else if ((this->error_region_prob[4] < vn) && (vn < this->error_region_prob[3])) {
                        soft1_pattern |= UNIT << bit_idx;
                    } else if ((this->error_region_prob[6] < vn) && (vn < this->error_region_prob[5])) {
                        soft0_pattern |= UNIT << bit_idx;
                    }
                }
            }

            cw.hard.cols[col_idx][word_idx] ^= error_pattern;
            cw.soft0.cols[col_idx][word_idx] ^= soft0_pattern;
            cw.soft1.cols[col_idx][word_idx] ^= soft1_pattern;
        }
    }

    return;
}

void error_injector::inject_fixed_count_errors(decoder_input_cw &cw) {
    const uint64_t UNIT = 1;
    uint64_t total_bits = (uint64_t)codeword_size * 64 * codeword::WORD_COUNT;
    uint64_t active_bits = (uint64_t)(bytes_of_userdata + bytes_of_parity) * 8;
    uint64_t num_errors_to_inject = static_cast<uint64_t>(std::round(rber * active_bits));

    if (num_errors_to_inject > total_bits)
        num_errors_to_inject = total_bits;

    std::vector<uint64_t> bit_indices(total_bits);
    std::iota(bit_indices.begin(), bit_indices.end(), 0);
    std::shuffle(bit_indices.begin(), bit_indices.end(), gen);

    uint64_t errors_injected = 0;
    for (uint64_t i = 0; (i < bit_indices.size()) && (errors_injected < num_errors_to_inject); ++i) {
        uint64_t linear_bit_pos = bit_indices[i];
        int col_idx = linear_bit_pos / (64 * codeword::WORD_COUNT);
        int word_idx = (linear_bit_pos % (64 * codeword::WORD_COUNT)) / 64;
        int bit_in_word_pos = linear_bit_pos % 64;
        int k = word_idx * 64 + bit_in_word_pos;
        bool blank_bit = !bit_active(col_idx, k);

        if (!blank_bit) {
            cw.hard.cols[col_idx][word_idx] ^= (UNIT << bit_in_word_pos);
            errors_injected++;
        }
    }

    return;
}

void error_injector::set_error_region_probs_tbench(double rber) {
    double target_cdf;
    double target_rber;
    double sigma;
    double mu = -1.0;

    this->error_region_prob[0] = rber;
    this->thresholds[0] = 0.0;

    if (strobes == 1) {
        this->thresholds[1] = 0.0;
        this->thresholds[2] = 0.0;
        this->thresholds[3] = 0.0;
        this->thresholds[4] = 0.0;
        this->thresholds[5] = 0.0;
        this->thresholds[6] = 0.0;
    } else {
        if (strobes == 3) {
            target_rber = this->target_distribution_rber;
            target_cdf = 1.0 - target_rber;
            sigma = find_sigma(this->thresholds[0], mu, target_cdf);
            target_rber *= 0.10;
            target_cdf = 1.0 - target_rber;
            this->thresholds[2] = find_x(mu, sigma, target_cdf);
            this->thresholds[4] = this->thresholds[2];
            this->thresholds[6] = this->thresholds[4];
        }

        if (strobes == 5) {
            target_rber = this->target_distribution_rber;
            target_cdf = 1.0 - target_rber;
            sigma = find_sigma(this->thresholds[0], mu, target_cdf);
            target_rber *= 0.05;
            target_cdf = 1.0 - target_rber;
            this->thresholds[4] = find_x(mu, sigma, target_cdf);
            this->thresholds[2] = 0.5 * this->thresholds[4];
            this->thresholds[6] = this->thresholds[4];
        }

        if (strobes == 7) {
            target_rber = this->target_distribution_rber;
            target_cdf = 1.0 - target_rber;
            sigma = find_sigma(this->thresholds[0], mu, target_cdf);
            target_rber *= 0.03;
            target_cdf = 1.0 - target_rber;
            this->thresholds[6] = find_x(mu, sigma, target_cdf);
            this->thresholds[2] = 0.333333 * this->thresholds[6];
            this->thresholds[4] = 0.666667 * this->thresholds[6];
        }

        this->thresholds[1] = -this->thresholds[2];
        this->thresholds[3] = -this->thresholds[4];
        this->thresholds[5] = -this->thresholds[6];
    }

    target_cdf = 1.0 - rber;
    sigma = find_sigma(this->thresholds[0], mu, target_cdf);

    this->error_region_prob[0] = rber;
    this->error_region_prob[1] = 1.0 - cdf(this->thresholds[1], mu, sigma);
    this->error_region_prob[2] = 1.0 - cdf(this->thresholds[2], mu, sigma);
    this->error_region_prob[3] = 1.0 - cdf(this->thresholds[3], mu, sigma);
    this->error_region_prob[4] = 1.0 - cdf(this->thresholds[4], mu, sigma);
    this->error_region_prob[5] = 1.0 - cdf(this->thresholds[5], mu, sigma);
    this->error_region_prob[6] = 1.0 - cdf(this->thresholds[6], mu, sigma);

    return;
}

void error_injector::set_error_region_probs_cmodel(double rber) {
    auto normalCDF = [&](double value) -> double { return 0.5 * erfc(-value * M_SQRT1_2); };

    auto normalInverseCDF = [&](double value) -> double {
        double x = 1.0;
        double y = 0.0;
        bool finished = 0;

        for (int i = 0; i < 10000; i++) {
            y = normalCDF(x);
            x = x + (0.25 * (value - y));
            finished = 1;
        }

        (void)finished;
        return x;
    };

    double delta[3] = {0, 0, 0};

    switch (strobes) {
    case 3:
        delta[0] = 0.350;
        break;
    case 5:
        delta[0] = 0.250;
        delta[1] = 0.500;
        break;
    case 7:
        delta[0] = 0.150;
        delta[1] = 0.300;
        delta[2] = 0.500;
        break;
    default:
        delta[0] = 0;
        delta[1] = 0;
        delta[2] = 0;
        break;
    }

    this->thresholds[0] = normalInverseCDF(1.0 - rber);
    this->thresholds[1] = this->thresholds[0] * (1.0 - delta[0]);
    this->thresholds[2] = this->thresholds[0] * (1.0 + delta[0]);
    this->thresholds[3] = this->thresholds[0] * (1.0 - delta[1]);
    this->thresholds[4] = this->thresholds[0] * (1.0 + delta[1]);
    this->thresholds[5] = this->thresholds[0] * (1.0 - delta[2]);
    this->thresholds[6] = this->thresholds[0] * (1.0 + delta[2]);

    this->error_region_prob[0] = rber;
    this->error_region_prob[1] = 1.0 - normalCDF(this->thresholds[1]);
    this->error_region_prob[2] = 1.0 - normalCDF(this->thresholds[2]);
    this->error_region_prob[3] = 1.0 - normalCDF(this->thresholds[3]);
    this->error_region_prob[4] = 1.0 - normalCDF(this->thresholds[4]);
    this->error_region_prob[5] = 1.0 - normalCDF(this->thresholds[5]);
    this->error_region_prob[6] = 1.0 - normalCDF(this->thresholds[6]);

    return;
}

double error_injector::cdf(double x, double mu, double sigma) {
    return 0.5 * (1.0 + std::erf((x - mu) / (sigma * std::sqrt(2.0))));
}

double error_injector::find_x(double mu, double sigma, double target_cdf) {
    int i;
    double x = 0.0;
    double cdf_x = 0.0;
    double error;

    i = 0;
    error = 1;

    while ((std::abs(error) >= 0.00000001) && (i < 10000)) {
        cdf_x = cdf(x, mu, sigma);
        error = target_cdf - cdf_x;
        x += (0.1 * error);
        i++;
    }

    return x;
}

double error_injector::find_sigma(double x, double mu, double target_cdf) {
    int i;
    double sigma;
    double cdf_x = 0.0;
    double error;

    sigma = 1.0;
    i = 0;
    error = 1;

    while ((std::abs(error) >= 0.00000001) && (i < 10000)) {
        cdf_x = cdf(x, mu, sigma);
        error = target_cdf - cdf_x;
        sigma *= (1.0 - (1.0 * error));
        i++;
    }

    return sigma;
}
