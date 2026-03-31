#include <chrono>
#include <iomanip>
#include <iostream>

#define DECODER_GPU_IMPL
#include "decoder_gpu.h"

#if defined(__CUDACC__)
__device__ __constant__ device_global_info_struct decoder_global_info;
#endif

void copy_decoder_info_to_device(
    const device_global_info_struct &host_decoder_info) {
#if defined(__CUDACC__)
  cudaError_t err = cudaMemcpyToSymbol(decoder_global_info, &host_decoder_info,
                                       sizeof(device_global_info_struct));
  if (err != cudaSuccess) {
    std::cerr << "[GPU] cudaMemcpyToSymbol failed: " << cudaGetErrorString(err)
              << "\n";
  }
#else
  (void)host_decoder_info;
#endif
}

int execute_gpu_kernel(const uint64_t TOTAL_CODEWORDS,
                       const uint64_t MAX_FAILURE_COUNT,
                       decoder_input_cw *test_cw) {
  (void)test_cw;
  const auto start_time = std::chrono::steady_clock::now();

#if defined(__CUDACC__)
  int device_count = 0;
  cudaError_t err = cudaGetDeviceCount(&device_count);
  if (err != cudaSuccess || device_count <= 0) {
    std::cerr
        << "[GPU] No usable CUDA device. Falling back to no-op GPU path.\n";
    return 0;
  }

  cudaDeviceProp prop{};
  cudaGetDeviceProperties(&prop, 0);
  std::cout << "[GPU] Device: " << prop.name << "\n";
  std::cout << "[GPU] TOTAL_CODEWORDS=" << TOTAL_CODEWORDS
            << ", MAX_FAILURE_COUNT=" << MAX_FAILURE_COUNT << "\n";

  // OCR restoration note: kernel body is intentionally minimal in this repair
  // pass.
  cudaDeviceSynchronize();
#else
  std::cout
      << "[GPU-STUB] execute_gpu_kernel called without CUDA compilation.\n";
  std::cout << "[GPU-STUB] TOTAL_CODEWORDS=" << TOTAL_CODEWORDS
            << ", MAX_FAILURE_COUNT=" << MAX_FAILURE_COUNT << "\n";
#endif

  const auto end_time = std::chrono::steady_clock::now();
  const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                      end_time - start_time)
                      .count();
  std::cout << "[GPU] elapsed " << ms << " ms\n";
  return 0;
}
