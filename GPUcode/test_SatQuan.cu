#include <cstdio>

#if __has_include(<cuda_fp16.h>)
#include <cuda_fp16.h>
#include <cuda_runtime.h>
#define SATQUAN_HAS_CUDA 1
#else
#define SATQUAN_HAS_CUDA 0
using half = float;
static inline float half2float(half v) { return v; }
static inline half float2half(float v) { return v; }
#endif

#define N 16

#if SATQUAN_HAS_CUDA
__device__ inline half Sat_Quan(half x, half Max_V, half Min_V, int m, int k) {
  float xf = half2float(x);
  xf = fminf(half2float(Max_V), fmaxf(half2float(Min_V), xf));
  int scale = 1 << k;
  int q = __float2int_rn(xf * scale);
  int mask = (1 << m) - 1;
  q &= mask;
  q = (q ^ (1 << (m - 1))) - (1 << (m - 1));
  return __float2half((float)q / scale);
}

__global__ void test_kernel(half *in, half *out, half max_v, half min_v, int m,
                            int k) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < N)
    out[i] = Sat_Quan(in[i], max_v, min_v, m, k);
}
#endif

int main() {
  half h_in[N];
  half h_out[N];
  for (int i = 0; i < N; i++) {
    float v = -1.0f + i * (2.0f / (N - 1));
    h_in[i] = float2half(v);
    h_out[i] = 0;
  }

#if SATQUAN_HAS_CUDA
  half *d_in = nullptr;
  half *d_out = nullptr;
  cudaMalloc(&d_in, N * sizeof(half));
  cudaMalloc(&d_out, N * sizeof(half));
  cudaMemcpy(d_in, h_in, N * sizeof(half), cudaMemcpyHostToDevice);

  int m = 6;
  int k = 3;
  half max_v = float2half(1.0f);
  half min_v = float2half(-1.0f);
  test_kernel<<<1, N>>>(d_in, d_out, max_v, min_v, m, k);
  cudaDeviceSynchronize();
  cudaMemcpy(h_out, d_out, N * sizeof(half), cudaMemcpyDeviceToHost);
  cudaFree(d_in);
  cudaFree(d_out);
#else
  for (int i = 0; i < N; i++)
    h_out[i] = h_in[i];
#endif

  printf("Input\tOutput\n");
  for (int i = 0; i < N; i++) {
    printf("%f\t%f\n", half2float(h_in[i]), half2float(h_out[i]));
  }
  return 0;
}
