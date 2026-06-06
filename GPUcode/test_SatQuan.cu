#include <cuda_fp16.h>
#include <stdio.h>

#define N 16

__device__ inline __half Sat_Quan(__half x, __half Max_V, __half Min_V, int m, int k) {
    float xf = __half2float(x);
    xf = fminf(__half2float(Max_V), fmaxf(__half2float(Min_V), xf));
    int scale = 1 << k;
    int q = __float2int_rn(xf * scale);
    int mask = (1 << m) - 1;
    q &= mask;
    q = (q ^ (1 << (m - 1))) - (1 << (m - 1));
    return __float2half((float)q / scale);
}

__global__ void test_kernel(__half *in, __half *out, __half max_v, __half min_v, int m, int k) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < N)
        out[i] = Sat_Quan(in[i], max_v, min_v, m, k);
}

int main() {
    __half h_in[N];
    __half h_out[N];
    __half *d_in;
    __half *d_out;

    for (int i = 0; i < N; i++) {
        float v = -4.0f + i * 0.5f;
        h_in[i] = __float2half(v);
    }

    cudaMalloc(&d_in, N * sizeof(__half));
    cudaMalloc(&d_out, N * sizeof(__half));
    cudaMemcpy(d_in, h_in, N * sizeof(__half), cudaMemcpyHostToDevice);

    int m = 6;
    int k = 3;
    __half max_v = __float2half(3.5f);
    __half min_v = __float2half(-3.5f);

    test_kernel<<<1, N>>>(d_in, d_out, max_v, min_v, m, k);
    cudaDeviceSynchronize();
    cudaMemcpy(h_out, d_out, N * sizeof(__half), cudaMemcpyDeviceToHost);

    printf("Input\tOutput\n");
    for (int i = 0; i < N; i++)
        printf("%f\t%f\n", __half2float(h_in[i]), __half2float(h_out[i]));

    cudaFree(d_in);
    cudaFree(d_out);

    return 0;
}
