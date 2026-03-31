#pragma once

#include <cstdint>
#include <cstring>
#include <iostream>

#include "decoder.h"

#if __has_include(<cuda_runtime.h>)
#include <cuda_runtime.h>
#define LDPC_GPU_HAS_CUDA_RUNTIME 1
#else
#define LDPC_GPU_HAS_CUDA_RUNTIME 0
#endif

#if __has_include(<cuda_fp16.h>)
#include <cuda_fp16.h>
#define LDPC_GPU_HAS_HALF 1
#else
#define LDPC_GPU_HAS_HALF 0
#endif

struct row_info {
  uint16_t actual_idx = 0;
  uint16_t mask_weight = 0;
  uint16_t mask_offset = 0;
  uint16_t element = 0;
};

struct device_global_info_struct {
  // H matrix info
  int16_t cols = 0;
  int16_t rows = 0;
  int16_t bits = LDPC_P;
  int16_t e_pre[LDPC_M][LDPC_N] = {{0}};

  // Decoder input params
  uint16_t nand_strobes = 0;
  uint16_t soft_bits = 0;
  uint16_t iteration_limit = 0;

  // Decoder params
  uint8_t VN_BITS = 8;
  bool post_process_en = true;
  uint16_t syndrome_weight_thr_qc = 16;
  uint16_t syndrome_weight_thr_post = 48;
  uint16_t likelihood_thr = 128;
  uint16_t post_ratio = 12;
  uint16_t likelihood_init_coef_all[8][4] = {{0}};
  uint16_t likelihood_init_fraction[3] = {0, 6, 6};

  // Error injection params
  float error_region_prob[7] = {0.0f};
  float vref[7] = {0.0f};
  uint8_t max = 0;
  uint8_t flip_thr = 0;
  bool prng_verilog_mode = false;

  // H-matrix related data
  uint16_t element[LDPC_M][LDPC_N] = {{0}};
  int8_t fade[LDPC_M][LDPC_N] = {{0}};
  int8_t occupied[LDPC_M][LDPC_N] = {{0}};
  uint16_t mask[LDPC_N][32] = {{0}}; // 32x16 = 512 bits

  int8_t extra_bytes_of_parity = 0;
  int8_t extra_bytes_of_userdata = 0;
  int32_t extra_bits_of_parity = 0;

  // finite params
  int finite_mode = 0;
  int finite_q_num = 0;
  int finite_c_num = 0;
  int finite_f_num = 0;
  float finite_q_max = 0.0f;
  float finite_q_min = 0.0f;
  float finite_r_max = 0.0f;
  float finite_r_min = 0.0f;
  float finite_c_max = 0.0f;
  float finite_c_min = 0.0f;
  float alpha = 0.0f;

  float llr_table[8] = {0.0f};
  float awgn_sigma = 0.0f;
};

#if defined(__CUDACC__) || defined(DECODER_GPU_IMPL)
void copy_decoder_info_to_device(
    const device_global_info_struct &host_decoder_info);
int execute_gpu_kernel(const uint64_t TOTAL_CODEWORDS,
                       const uint64_t MAX_FAILURE_COUNT,
                       decoder_input_cw *test_cw = nullptr);
#else
static inline void copy_decoder_info_to_device(
    const device_global_info_struct &host_decoder_info) {
  (void)host_decoder_info;
}

static inline int execute_gpu_kernel(const uint64_t TOTAL_CODEWORDS,
                                     const uint64_t MAX_FAILURE_COUNT,
                                     decoder_input_cw *test_cw = nullptr) {
  (void)test_cw;
  std::cout << "[GPU-STUB] CUDA runtime unavailable in this build.\n";
  std::cout << "[GPU-STUB] Requested TOTAL_CODEWORDS=" << TOTAL_CODEWORDS
            << ", MAX_FAILURE_COUNT=" << MAX_FAILURE_COUNT << "\n";
  return 0;
}
#endif
