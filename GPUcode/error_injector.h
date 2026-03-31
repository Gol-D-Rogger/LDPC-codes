#pragma once

#include <cstdint>
#include <functional>
#include <random>

#include "codeword.h"

struct error_injector {
  static constexpr double DEFAULT_TARGET_DIST_RBER = 5.0;
  static constexpr bool RANDOM_SEED = false;
  static constexpr int MAX_THRESHOLDS = 7;

  int insertion_mode = 0;
  int prob_region_mode = 0;
  int strobes = 0;
  double target_distribution_rber = DEFAULT_TARGET_DIST_RBER;
  int num_soft_bits = 0;
  double rber = 0;
  int parity_col_idx = 0;
  int final_parity_col_bits = 0;
  int final_udata_col_bits = 0;
  int codeword_size = 0;

  double thresholds[MAX_THRESHOLDS] = {0.0};
  double error_region_prob[MAX_THRESHOLDS] = {0.0};

  uint64_t seed = 784225534;
  std::random_device rd;
  std::mt19937 gen;

  error_injector();
  error_injector(int insertion_mode, int strobes, double rber, int soft_bits,
                 int parity_col_idx, int final_udata_col_bits,
                 int final_parity_col_bits, int codeword_size,
                 double target_distribution_rber = DEFAULT_TARGET_DIST_RBER);

  void set_error_region_probs_tbench(double rber);
  void set_error_region_probs_cmodel(double rber);
  void inject_normal_distribution_errors(decoder_input_cw &cw);
  void inject_fixed_count_errors(decoder_input_cw &cw);

  std::function<void(decoder_input_cw &)> inject_errors;

  double find_x(double mu, double sigma, double target_cdf);
  double find_sigma(double x, double mu, double target_cdf);
  double cdf(double x, double mu, double sigma);
};
