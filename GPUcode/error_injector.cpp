#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "error_injector.h"

namespace {
double clamp01(double v) {
  if (v < 0.0)
    return 0.0;
  if (v > 1.0)
    return 1.0;
  return v;
}
} // namespace

error_injector::error_injector() {
  seed = RANDOM_SEED ? rd() : seed;
  gen.seed(seed);
  inject_errors = std::bind(&error_injector::inject_normal_distribution_errors,
                            this, std::placeholders::_1);
}

error_injector::error_injector(int insertion_mode, int strobes, double rber,
                               int soft_bits, int parity_col_idx,
                               int final_udata_col_bits,
                               int final_parity_col_bits, int codeword_size,
                               double target_distribution_rber) {
  this->insertion_mode = insertion_mode;
  this->strobes = strobes;
  this->rber = rber;
  this->parity_col_idx = parity_col_idx;
  this->final_udata_col_bits = final_udata_col_bits;
  this->final_parity_col_bits = final_parity_col_bits;
  this->codeword_size = codeword_size;
  this->target_distribution_rber = target_distribution_rber;
  this->num_soft_bits = soft_bits;

  seed = RANDOM_SEED ? rd() : seed;
  gen.seed(seed);

  if (insertion_mode == 0) {
    inject_errors =
        std::bind(&error_injector::inject_normal_distribution_errors, this,
                  std::placeholders::_1);
  } else {
    inject_errors = std::bind(&error_injector::inject_fixed_count_errors, this,
                              std::placeholders::_1);
  }

  if (prob_region_mode == 0)
    set_error_region_probs_tbench(rber);
  else
    set_error_region_probs_cmodel(rber);
}

void error_injector::set_error_region_probs_tbench(double rber) {
  const double safe_rber = clamp01(rber);
  error_region_prob[0] = safe_rber;

  double spread1 = 0.5;
  double spread2 = 1.0;
  double spread3 = 1.5;
  if (strobes <= 1) {
    spread1 = 0.0;
    spread2 = 0.0;
    spread3 = 0.0;
  } else if (strobes == 3) {
    spread2 = spread1;
    spread3 = spread1;
  } else if (strobes == 5) {
    spread3 = spread2;
  }

  thresholds[0] = 0.0;
  thresholds[1] = -spread1;
  thresholds[2] = +spread1;
  thresholds[3] = -spread2;
  thresholds[4] = +spread2;
  thresholds[5] = -spread3;
  thresholds[6] = +spread3;

  const double mu = 0.0;
  const double sigma = 1.0;
  error_region_prob[1] = 1.0 - cdf(thresholds[1], mu, sigma);
  error_region_prob[2] = 1.0 - cdf(thresholds[2], mu, sigma);
  error_region_prob[3] = 1.0 - cdf(thresholds[3], mu, sigma);
  error_region_prob[4] = 1.0 - cdf(thresholds[4], mu, sigma);
  error_region_prob[5] = 1.0 - cdf(thresholds[5], mu, sigma);
  error_region_prob[6] = 1.0 - cdf(thresholds[6], mu, sigma);

  for (int i = 0; i < MAX_THRESHOLDS; ++i)
    error_region_prob[i] = clamp01(error_region_prob[i]);
}

void error_injector::set_error_region_probs_cmodel(double rber) {
  // Keep a conservative fallback aligned with current OCR intent.
  set_error_region_probs_tbench(rber);
}

void error_injector::inject_normal_distribution_errors(decoder_input_cw &cw) {
  const uint64_t UNIT = 1ULL;
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  const int col_limit =
      std::max(0, std::min(codeword_size, codeword::MAX_COL_COUNT));

  for (int col_idx = 0; col_idx < col_limit; ++col_idx) {
    for (int word_idx = 0; word_idx < codeword::WORD_COUNT; ++word_idx) {
      uint64_t error_pattern = 0;
      uint64_t soft0_pattern = 0;
      uint64_t soft1_pattern = 0;

      for (unsigned bit_idx = 0; bit_idx < 64; ++bit_idx) {
        const int k = word_idx * 64 + static_cast<int>(bit_idx);
        bool blank_bit = false;
        blank_bit = blank_bit || ((final_udata_col_bits > 0) &&
                                  (col_idx == parity_col_idx - 1) &&
                                  (k >= final_udata_col_bits));
        blank_bit = blank_bit || ((final_parity_col_bits > 0) &&
                                  (col_idx == parity_col_idx) &&
                                  (k >= final_parity_col_bits));
        if (blank_bit)
          continue;

        const double vn = dist(gen);
        if (vn < rber)
          error_pattern |= (UNIT << bit_idx);

        if (num_soft_bits == 1) {
          if ((error_region_prob[2] < vn) && (vn < error_region_prob[1])) {
            soft0_pattern |= (UNIT << bit_idx);
            soft1_pattern |= (UNIT << bit_idx);
          }
        } else if (num_soft_bits >= 2) {
          if ((error_region_prob[2] < vn) && (vn < error_region_prob[1])) {
            soft0_pattern |= (UNIT << bit_idx);
            soft1_pattern |= (UNIT << bit_idx);
          } else if ((error_region_prob[4] < vn) &&
                     (vn < error_region_prob[3])) {
            soft1_pattern |= (UNIT << bit_idx);
          } else if ((error_region_prob[6] < vn) &&
                     (vn < error_region_prob[5])) {
            soft0_pattern |= (UNIT << bit_idx);
          }
        }
      }

      cw.hard.cols[col_idx][word_idx] ^= error_pattern;
      cw.soft0.cols[col_idx][word_idx] ^= soft0_pattern;
      cw.soft1.cols[col_idx][word_idx] ^= soft1_pattern;
    }
  }
}

void error_injector::inject_fixed_count_errors(decoder_input_cw &cw) {
  const uint64_t UNIT = 1ULL;
  const int col_limit =
      std::max(0, std::min(codeword_size, codeword::MAX_COL_COUNT));
  const uint64_t total_bits =
      static_cast<uint64_t>(col_limit) * 64ULL * codeword::WORD_COUNT;
  uint64_t num_errors_to_inject = static_cast<uint64_t>(
      std::llround(rber * static_cast<double>(total_bits)));
  if (num_errors_to_inject > total_bits)
    num_errors_to_inject = total_bits;

  std::vector<uint64_t> bit_indices(total_bits);
  std::iota(bit_indices.begin(), bit_indices.end(), 0);
  std::shuffle(bit_indices.begin(), bit_indices.end(), gen);

  uint64_t errors_injected = 0;
  for (uint64_t i = 0;
       (i < bit_indices.size()) && (errors_injected < num_errors_to_inject);
       ++i) {
    const uint64_t linear_bit_pos = bit_indices[i];
    const int col_idx =
        static_cast<int>(linear_bit_pos / (64ULL * codeword::WORD_COUNT));
    const int word_idx = static_cast<int>(
        (linear_bit_pos % (64ULL * codeword::WORD_COUNT)) / 64ULL);
    const int bit_in_word_pos = static_cast<int>(linear_bit_pos % 64ULL);

    const int k = word_idx * 64 + bit_in_word_pos;
    bool blank_bit = false;
    blank_bit = blank_bit || ((final_udata_col_bits > 0) &&
                              (col_idx == parity_col_idx - 1) &&
                              (k >= final_udata_col_bits));
    blank_bit = blank_bit ||
                ((final_parity_col_bits > 0) && (col_idx == parity_col_idx) &&
                 (k >= final_parity_col_bits));
    if (blank_bit)
      continue;

    cw.hard.cols[col_idx][word_idx] ^= (UNIT << bit_in_word_pos);
    errors_injected++;
  }
}

double error_injector::cdf(double x, double mu, double sigma) {
  const double z = (x - mu) / (sigma * std::sqrt(2.0));
  return 0.5 * (1.0 + std::erf(z));
}

double error_injector::find_x(double mu, double sigma, double target_cdf) {
  double x = mu;
  for (int i = 0; i < 10000; ++i) {
    const double cdf_x = cdf(x, mu, sigma);
    const double error = target_cdf - cdf_x;
    if (std::abs(error) < 1e-6)
      break;
    x += 0.1 * error * sigma;
  }
  return x;
}

double error_injector::find_sigma(double x, double mu, double target_cdf) {
  double sigma = 1.0;
  for (int i = 0; i < 10000; ++i) {
    if (sigma < 1e-9)
      sigma = 1e-9;
    const double cdf_x = cdf(x, mu, sigma);
    const double error = target_cdf - cdf_x;
    if (std::abs(error) < 1e-6)
      break;
    sigma *= (1.0 - 0.1 * error);
  }
  return sigma;
}
