#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdio.h>

#include "decoder.h"
#include "decoder_gpu.h"

#define WARP_SIZE 32
#define NUM_WARPS 1
#define P_SIZE 512
#define P_SIZEm1 511
#define OUTPUT_BUFFER_SIZE 512
#define FINITE_MODE

//#define FAST_QUANT
//#define COMPARE_CPU

#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

__constant__ device_global_info_struct decoder_global_info;
__device__ curandState *d_states;

__global__ void ldpc_decode_kernel(decoder_input_cw *input_cw, ldpc_decoder_output *decoder_output, uint32_t *iter_info, decoder_workspace *cw_workspace, int bypass_errinj);
__global__ void init_rng(unsigned long seed);

__device__ void inject_normal_errors(decoder_input_cw *input_cw);
__device__ void inject_normal_errors2(decoder_input_cw *input_cw);
__device__ void config_check_nodes(const codeword *input_cw, check_nodes_gpu *cn);
__device__ void compute_likelihood_levels(int syndrome_weight, OptimizedSharedMemory *cw_shared_mem);
__device__ inline uint32_t check_node_weight(const check_nodes_gpu *cn);
__device__ inline uint32_t warp_reduce_sum(uint32_t val);
__device__ inline uint16_t mask_range_512_cuda(int offset, int count);
__device__ inline uint16_t mask_range_512_cuda_rot(int offset, int count);
__device__ inline uint16_t rotate_right_512(const uint16_t word, unsigned r);
__device__ inline uint16_t update_vn_planar(const uint16_t weights, uint32_t *vn_likelihood, const OptimizedSharedMemory *cw_shared_mem, const uint8_t post_triggers, const int bit_offset);
//__device__ inline uint16_t rotate_right_512_deng(uint16_t word, unsigned r);
__device__ inline void update_hd(int col, int layer, int layer_pre, short &hd_updated, decoder_workspace *cw_workspace, OptimizedSharedMemory *cw_shared_mem, int &shift_val1);
__device__ inline void update_node_next(int cir_cnt, int layer, int col, OptimizedSharedMemory *cw_shared_mem, decoder_workspace *cw_workspace, int shift_val1);
__device__ inline void update_node(int layer_pre, int col, decoder_workspace *cw_workspace, OptimizedSharedMemory *cw_shared_mem, int layer);
__device__ inline uint32_t calculate_crc32(const codeword *input_cw);
__device__ inline void print_512float(__half *a);
__device__ inline void print_512cn_msg(cn_msg *a);

__device__ inline void print_cn_msg(cn_msg cn) {
    printf("%f %f %d %d\n", __half2float(cn.min1_val), __half2float(cn.min2_val), cn.min1_pos, cn.sign_tot);
}

#define MEMBER_SIZE(type, member) sizeof(((type *)0)->member)

template <typename Kernel> int detail_kernel_occupancy(Kernel kernel, const char *kernel_name) {
    std::cout << "\n============== KERNEL OCCUPANCY DETAILS ==============\n";
    std::cout << ":: Kernel name: " << kernel_name << "\n";

    cudaFuncAttributes attr;
    int maxActiveBlocks;
    CUDA_CHECK(cudaFuncGetAttributes(&attr, kernel));
    CUDA_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&maxActiveBlocks, kernel, NUM_WARPS * WARP_SIZE, 0));

    int warps_per_block = (NUM_WARPS * WARP_SIZE + 31) / 32;
    int active_warps_per_SM = maxActiveBlocks * warps_per_block;

    std::cout << " - [" << kernel_name << "] Registers per thread: " << attr.numRegs << "\n";
    std::cout << " - [" << kernel_name << "] Max threads per block: " << attr.maxThreadsPerBlock << "\n";
    std::cout << " - [" << kernel_name << "] USED shared mem per block: " << attr.sharedSizeBytes << "\n";
    std::cout << " - [" << kernel_name << "] ACTIVE warps per block: " << warps_per_block << "\n";
    std::cout << " - [" << kernel_name << "] ACTIVE blocks per SM: " << maxActiveBlocks << "\n";
    std::cout << " - [" << kernel_name << "] ACTIVE warps per SM: " << active_warps_per_SM << "\n";

    int device = 0;
    CUDA_CHECK(cudaSetDevice(device));

    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, device));

    int warp_size = prop.warpSize;
    int max_threads_per_SM = prop.maxThreadsPerMultiProcessor;
    int max_warps_per_SM = max_threads_per_SM / warp_size;

    std::cout << "\n================= CUDA DEVICE INFO =================\n";
    std::cout << " - Device name: " << prop.name << "\n";
    std::cout << " - [" << prop.name << "] Total SMs: " << prop.multiProcessorCount << "\n";
    std::cout << " - [" << prop.name << "] Max threads per SM: " << max_threads_per_SM << "\n";
    std::cout << " - [" << prop.name << "] Max warps per SM: " << max_warps_per_SM << "\n";
    std::cout << " - [" << prop.name << "] Max blocks per SM: " << prop.maxBlocksPerMultiProcessor << "\n";
    std::cout << " - [" << prop.name << "] Registers per SM: " << prop.regsPerMultiprocessor << "\n";
    std::cout << " - [" << prop.name << "] Shared memory per SM: " << prop.sharedMemPerMultiprocessor << " bytes\n";
    std::cout << " - [" << prop.name << "] Shared memory per Block: " << prop.sharedMemPerBlock << " bytes\n";

    return prop.multiProcessorCount;
}

__host__ void copy_decoder_info_to_device(const device_global_info_struct &host_decoder_info) {
    CUDA_CHECK(cudaMemcpyToSymbol(decoder_global_info, &host_decoder_info, sizeof(device_global_info_struct)));
    logger::ITER_LIMIT = host_decoder_info.iteration_limit;
}

__host__ int execute_gpu_kernel(const uint64_t TOTAL_CODEWORDS, const uint64_t MAX_FAILURE_COUNT, decoder_input_cw *test_cw) {
    decoder_input_cw *input_cw = NULL;
    ldpc_decoder_output *decoder_output = NULL;
    decoder_workspace *cw_workspace = NULL;
    curandState *curand_states_ptr = NULL;
    uint32_t *iter_info = NULL;
    int bypass_errinj = 0;
    const int SM_COUNT = detail_kernel_occupancy(ldpc_decode_kernel, "ldpc_decode_kernel");

    const uint64_t MAX_MEM_CODEWORDS = SM_COUNT * 40;
    uint64_t BATCH_SIZE = std::min(TOTAL_CODEWORDS, MAX_MEM_CODEWORDS);
    if (BATCH_SIZE == 0)
        BATCH_SIZE = 1;

    const int NUM_BLOCKS = (BATCH_SIZE + NUM_WARPS - 1) / NUM_WARPS;
    const int ACTUAL_BATCH_SIZE = NUM_BLOCKS * NUM_WARPS;

    std::cout << " :: SM COUNT: " << SM_COUNT << ", Batch Size: " << ACTUAL_BATCH_SIZE << ", NUM_BLOCKS: " << NUM_BLOCKS << "\n";

    CUDA_CHECK(cudaMalloc((void **)&decoder_output, sizeof(ldpc_decoder_output) * OUTPUT_BUFFER_SIZE));
    CUDA_CHECK(cudaMalloc((void **)&cw_workspace, sizeof(decoder_workspace) * ACTUAL_BATCH_SIZE));
    CUDA_CHECK(cudaMalloc((void **)&input_cw, sizeof(decoder_input_cw) * ACTUAL_BATCH_SIZE));
    CUDA_CHECK(cudaMalloc((void **)&iter_info, sizeof(uint32_t) * logger::MAX_ITER));

    ldpc_decoder_output *host_decoder_output = (ldpc_decoder_output *)calloc(OUTPUT_BUFFER_SIZE, sizeof(*host_decoder_output));

    if (test_cw != NULL) {
        bypass_errinj = 1;
        for (int i = 0; i < ACTUAL_BATCH_SIZE; i++)
            CUDA_CHECK(cudaMemcpy(&input_cw[i], test_cw, sizeof(decoder_input_cw), cudaMemcpyHostToDevice));
    } else {
        std::cout << "\n======================= Initializing CUDA RNG =======================\n";
        CUDA_CHECK(cudaMalloc(&curand_states_ptr, ACTUAL_BATCH_SIZE * WARP_SIZE * sizeof(curandState)));
        CUDA_CHECK(cudaMemcpyToSymbol(d_states, &curand_states_ptr, sizeof(curandState *)));

        unsigned long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        init_rng<<<NUM_BLOCKS, WARP_SIZE * NUM_WARPS>>>(seed);
        CUDA_CHECK(cudaDeviceSynchronize());
        std::cout << " :: RNG initialized\n";
    }

    std::cout << "\n======================= Launching decode kernel =======================\n";

    const int TOTAL_BATCH_COUNT = (TOTAL_CODEWORDS + ACTUAL_BATCH_SIZE - 1) / ACTUAL_BATCH_SIZE;
    int BATCHES_EXECUTED = 0;
    auto start_time = std::chrono::steady_clock::now();

    decoder_output_acc final_result;
    memset(&final_result, 0, sizeof(final_result));

    for (BATCHES_EXECUTED = 0; (BATCHES_EXECUTED < TOTAL_BATCH_COUNT) && (final_result.failure < MAX_FAILURE_COUNT); BATCHES_EXECUTED++) {
        CUDA_CHECK(cudaMemset(cw_workspace, 0, sizeof(decoder_workspace) * ACTUAL_BATCH_SIZE));
        CUDA_CHECK(cudaMemset(decoder_output, 0, sizeof(ldpc_decoder_output) * OUTPUT_BUFFER_SIZE));

        ldpc_decode_kernel<<<NUM_BLOCKS, WARP_SIZE * NUM_WARPS>>>(input_cw, decoder_output, iter_info, cw_workspace, bypass_errinj);
        CUDA_CHECK(cudaDeviceSynchronize());
        CUDA_CHECK(cudaMemcpy(host_decoder_output, decoder_output, sizeof(ldpc_decoder_output) * OUTPUT_BUFFER_SIZE, cudaMemcpyDeviceToHost));

        for (int cw_id = 0; cw_id < OUTPUT_BUFFER_SIZE; cw_id++) {
            final_result.failure += host_decoder_output[cw_id].failure;
            final_result.total_errors += host_decoder_output[cw_id].total_errors;
            final_result.iterations += host_decoder_output[cw_id].iterations;
            final_result.clock_cycles += host_decoder_output[cw_id].clock_cycles;
            final_result.syndrome_weight_before += host_decoder_output[cw_id].syndrome_weight_before;
            final_result.syndrome_weight_after += host_decoder_output[cw_id].syndrome_weight_after;
            final_result.net_acc += host_decoder_output[cw_id].net_acc;
        }

        if (BATCHES_EXECUTED % 20 == 19) {
            std::cout << "======== " << (BATCHES_EXECUTED + 1) << "/" << TOTAL_BATCH_COUNT << " Batches Completed | ";
            std::cout << final_result.net_acc << " CWs | " << final_result.failure << " failures | ";
            if (final_result.net_acc > 0)
                std::cout << (double)final_result.iterations / final_result.net_acc << " avg. iters =======\n";
            else
                std::cout << "0 avg. iters ======\n";
        }

        if (BATCHES_EXECUTED % 400 == 0) {
            auto end_time = std::chrono::steady_clock::now();
            logger::print_accumulated_stats(final_result.net_acc, final_result);
            logger::log_elapsed_time(start_time, end_time);
        }
    }

    auto end_time = std::chrono::steady_clock::now();
    logger::log_elapsed_time(start_time, end_time);

    std::cout << "All execution done\n";
    logger::print_accumulated_stats(final_result.net_acc, final_result);
    if (final_result.net_acc > 0) {
        std::cout << "Overall Decoding Failure Rate: " << std::fixed << std::setprecision(10) << (double)final_result.failure / final_result.net_acc << "(" << final_result.failure << "/" << final_result.net_acc << ")\n";
    } else {
        std::cout << "Overall Decoding Failure Rate: 0 (0 / 0)\n";
    }

    uint32_t *host_iter_info = (uint32_t *)malloc(sizeof(uint32_t) * logger::MAX_ITER);
    CUDA_CHECK(cudaMemcpy(host_iter_info, iter_info, sizeof(uint32_t) * logger::MAX_ITER, cudaMemcpyDeviceToHost));
    logger::print_iter_stats(host_iter_info);

    CUDA_CHECK(cudaFree(input_cw));
    CUDA_CHECK(cudaFree(decoder_output));
    CUDA_CHECK(cudaFree(cw_workspace));
    if (curand_states_ptr)
        CUDA_CHECK(cudaFree(curand_states_ptr));

    free(host_decoder_output);
    free(host_iter_info);

    std::cout << "\n======================= Exiting GPU section =======================\n";
    return 0;
}

#ifdef FAST_QUANT
__device__ inline __half Sat_Quan(__half x, __half Max_V, __half Min_V, int m, int k) {
    float xf = __half2float(x);

    xf = fminf(__half2float(Max_V), fmaxf(__half2float(Min_V), xf));

    int scale = 1 << k;
    int q = __float2int_rn(xf * scale);

    int mask = (1 << m) - 1;
    int sign_bit = 1 << (m - 1);
    q &= mask;
    q = (q ^ sign_bit) - sign_bit;
    return __float2half((float)q / scale);
}
#else
__device__ inline int Quantize(__half x, int m, int k) {
    int t, t1;
    t = 1 << k;
    t = (int)(x * __int2half_rn(t) + (x > __half(0.0) ? __half(0.5) : __half(-0.5)));
    t1 = (1 << m) - 1;
    t1 = t1 & t;
    return t1;
}

__device__ inline void Tru2intS(int x, int m, int *dx) {
    int t1 = x;
    int t2 = -1;
    t1 = (x >> (m - 1));
    if (t1 == 0)
        (*dx) = x;
    else {
        t1 = (1 << m) - 1;
        t2 = t2 ^ t1; // ^ is 'xor' operation;
        *dx = t2 | x;
    }
}

__device__ inline __half Sat_Quan(__half x, __half Max_V, __half Min_V, int m, int k) {
    __half val;
    int q, p;

    val = x;
    val = val > Max_V ? Max_V : val;
    val = val < Min_V ? Min_V : val;

    q = Quantize(val, m, k);
    Tru2intS(q, m, &p);
    return p / pow(2, k);
}

#endif

__global__ void init_rng(unsigned long seed) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    curand_init(seed, idx, 0, &d_states[idx]);
}

__global__ void ldpc_decode_kernel(decoder_input_cw *input_cw, ldpc_decoder_output *decoder_output, uint32_t *iter_info, decoder_workspace *cw_workspace, int bypass_errinj) {
    const int warp_id = threadIdx.x / 32;
    const int lane_id = threadIdx.x & 31;
    const int cw_id = blockIdx.x * NUM_WARPS + warp_id;

    short ldec_early_term_en = 1;
    short hd_stable_cnt = 0;
    short synd_pass_cnt = 0;
    short cir_cnt = 0;
    short cw_fail = 0;
    short hd_init = 0;
    short hd_updated = 0;

    __shared__ OptimizedSharedMemory cw_shared_mem[NUM_WARPS];

    __half alpha = __float2half(decoder_global_info.alpha);
    bool finite_mode = decoder_global_info.finite_mode;

    if (bypass_errinj == 0) {
#ifdef USE_REGION_INJECT
        inject_normal_errors(&input_cw[cw_id]);
#else
        inject_normal_errors2(&input_cw[cw_id]);
#endif
    }

    if (lane_id == 0) {
        cw_shared_mem[warp_id].syndrome_weight = 0;
        cw_shared_mem[warp_id].finished = 0;
    }

    for (int idx = lane_id; idx < decoder_global_info.cols * 8; idx += WARP_SIZE) {
        const auto col_idx = idx / 8;
        const auto word_idx = idx % 8;
        const auto hard_bits = input_cw[cw_id].hard.cols[col_idx][word_idx];
        const auto soft_bits_0 = input_cw[cw_id].soft0.cols[col_idx][word_idx];
        const auto soft_bits_1 = input_cw[cw_id].soft1.cols[col_idx][word_idx];

        cw_workspace[cw_id].vn.dec_do_blk.cols[col_idx][word_idx] = hard_bits;

        for (int bit_idx = 0; bit_idx < 64; bit_idx++) {
            const bool hard = (hard_bits & (1UL << bit_idx));
            const bool sb0 = (soft_bits_0 & (1UL << bit_idx));
            const bool sb1 = (soft_bits_1 & (1UL << bit_idx));
            const auto llr_index = (hard << 2) | (sb1 << 1) | sb0;
            cw_workspace[cw_id].vn.cn_q_mem[col_idx][word_idx * 64 + bit_idx] = __float2half(decoder_global_info.llr_table[llr_index]);
        }
    }

    for (int i = lane_id; i < decoder_global_info.cols; i += WARP_SIZE)
        cw_shared_mem[warp_id].dec_init[i] = 1;

    __syncwarp();

    int iteration = 0;
    while (iteration < decoder_global_info.iteration_limit && !cw_shared_mem[warp_id].finished) {
        cir_cnt = 0;
        cw_fail = 1;

        for (int layer = 0; layer < decoder_global_info.rows && ((ldec_early_term_en == 0) || (cw_fail == 1)); layer++) {
            hd_init = 0;
            for (int i = 0; i < decoder_global_info.cols; i++) {
                if (cw_shared_mem[warp_id].dec_init[i] != 0) {
                    hd_init = 1;
                    break;
                }
            }

            for (int i = lane_id; i < P_SIZE; i += WARP_SIZE) {
                cw_workspace[cw_id].cn_c_updt_cur[i].min1_val = __float2half(10000.0f);
                cw_workspace[cw_id].cn_c_updt_cur[i].min2_val = __float2half(10000.0f);
                cw_workspace[cw_id].cn_c_updt_cur[i].min1_pos = 0;
                cw_workspace[cw_id].cn_c_updt_cur[i].sign_tot = 1;
            }

            hd_updated = 0;
            cw_shared_mem[warp_id].layer_synd[lane_id] = 0;
            __syncwarp();

            for (int col = 0; col < decoder_global_info.cols; col++) {
                if (((decoder_global_info.occupied[layer][col] == 0) && (decoder_global_info.fade[layer][col] == 0)) || ((decoder_global_info.fade[layer][col] == 1) && (decoder_global_info.extra_bytes_of_parity == 0)))
                    continue;

                int layer_pre = decoder_global_info.e_pre[layer][col];
                int shift_val1 = 0;

                update_node(layer_pre, col, cw_workspace, &cw_shared_mem[warp_id], layer);
                update_hd(col, layer, layer_pre, hd_updated, cw_workspace, &cw_shared_mem[warp_id], shift_val1);
                update_node_next(cir_cnt, layer, col, &cw_shared_mem[warp_id], cw_workspace, shift_val1);
                cir_cnt++;
            } // per circulant

// update C_MSG per layer
#ifdef USE_CN_MSG
            auto *dst = cw_workspace[cw_id].vn.cn_c_mem[layer];
#else
            auto *dst_min1_val = cw_workspace[cw_id].vn.cn_c_mem_min1_val[layer];
            auto *dst_min2_val = cw_workspace[cw_id].vn.cn_c_mem_min2_val[layer];
            auto *dst_min1_pos = cw_workspace[cw_id].vn.cn_c_mem_min1_pos[layer];
#endif

            auto *src = cw_workspace[cw_id].cn_c_updt_cur;
            for (int i = lane_id; i < P_SIZE; i += WARP_SIZE) {
                __half min1 = src[i].min1_val * alpha;
                __half min2 = src[i].min2_val * alpha;

#ifdef FINITE_MODE
                min1 = Sat_Quan(min1, decoder_global_info.finite_c_max, decoder_global_info.finite_c_min, decoder_global_info.finite_c_num, decoder_global_info.finite_f_num);
                min2 = Sat_Quan(min2, decoder_global_info.finite_c_max, decoder_global_info.finite_c_min, decoder_global_info.finite_c_num, decoder_global_info.finite_f_num);
#endif
#ifdef USE_CN_MSG
                cn_msg tmp_c;
                tmp_c.min1_val = min1;
                tmp_c.min2_val = min2;
                tmp_c.min1_pos = src[i].min1_pos;
                tmp_c.sign_tot = src[i].sign_tot;
                dst[i] = tmp_c;
#else
                dst_min1_val[i] = min1;
                dst_min2_val[i] = min2;
                dst_min1_pos[i] = src[i].min1_pos + ((src[i].sign_tot == 1) ? 128 : 0);
#endif
            }

            uint64_t *p = (uint64_t *)cw_shared_mem[warp_id].layer_synd;
            int layer_synd_wt = 0;
#pragma unroll
            for (int i = 0; i < 8; i++)
                layer_synd_wt += __popcll(p[i]);

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

            if ((synd_pass_cnt >= decoder_global_info.rows) && (hd_stable_cnt >= (decoder_global_info.rows - 1))) {
                cw_fail = 0;
                if (lane_id == 0)
                    cw_shared_mem[warp_id].finished = 1;
            }
        }

        iteration++;
        __syncwarp();
    }

    if (lane_id == 0) {
        atomicAdd(&decoder_output[cw_id % OUTPUT_BUFFER_SIZE].failure, (cw_fail != 0));
        atomicAdd(&decoder_output[cw_id % OUTPUT_BUFFER_SIZE].iterations, iteration);
        atomicAdd(&decoder_output[cw_id % OUTPUT_BUFFER_SIZE].net_acc, 1);
        atomicAdd(&iter_info[iteration], 1);
    }
    __syncwarp();

    return;
}

__device__ void config_check_nodes(const codeword *input_cw, check_nodes_gpu *cn) {
    const int lane_id = threadIdx.x & 31;

    for (int i = lane_id; i < check_nodes_gpu::MAX_ROW_COUNT * 8; i += WARP_SIZE) {
        const auto row_idx = i / 8;
        const auto word_idx = i % 8;
        cn->rows[row_idx][word_idx] = 0ULL;
    }

    __syncwarp();

    for (unsigned col_idx = 0; col_idx < decoder_global_info.cols; col_idx++) {
        const auto active_row_count = decoder_global_info.active_rows_in_col[col_idx];
        for (unsigned word_idx = 0; word_idx < codeword::WORD_COUNT; word_idx++) {
            const uint64_t word = input_cw->cols[col_idx][word_idx];
            if (word == 0)
                continue;

            for (unsigned bit_idx_in_word = 0; bit_idx_in_word < codeword::WORD_SIZE; bit_idx_in_word++) {
                const unsigned bit_idx = word_idx * codeword::WORD_SIZE + bit_idx_in_word;

                if (word & (codec_support::UNIT << bit_idx_in_word)) {
                    for (unsigned row_idx = lane_id; row_idx < active_row_count; row_idx += WARP_SIZE) {
                        const auto actual_row_idx = decoder_global_info.elem_info_trnsp[col_idx][row_idx].actual_idx;
                        const auto active_count = decoder_global_info.elem_info_trnsp[col_idx][row_idx].mask_weight;
                        const auto offset = decoder_global_info.elem_info_trnsp[col_idx][row_idx].mask_offset;
                        const auto element = decoder_global_info.elem_info_trnsp[col_idx][row_idx].element;

                        unsigned bit_loc = (bit_idx - offset) & 0x1FF;

                        if (bit_loc < active_count) {
                            const uint32_t rotated_bit_idx = (bit_idx - element) & 0x1FF;
                            const uint32_t rotated_word_idx = rotated_bit_idx / 64;
                            const uint32_t rotated_bit_idx_in_word = rotated_bit_idx % 64;

                            cn->rows[actual_row_idx][rotated_word_idx] ^= (codec_support::UNIT << rotated_bit_idx_in_word);
                        }
                    }
                }
            }
        }
    }

    __syncwarp();
    return;
}

__device__ uint32_t check_node_weight(const check_nodes_gpu *cn) {
    const int lane_id = threadIdx.x & 31;
    uint32_t local_weight = 0;
    const auto idx_limit = decoder_global_info.rows * check_nodes_gpu::WORD_COUNT;

    for (int idx = lane_id; idx < idx_limit; idx += WARP_SIZE) {
        const auto row_idx = idx / check_nodes_gpu::WORD_COUNT;
        const auto word_idx = idx % check_nodes_gpu::WORD_COUNT;
        local_weight += __popcll(cn->rows[row_idx][word_idx]);
    }

    return warp_reduce_sum(local_weight);
}

__device__ inline uint32_t warp_reduce_sum(uint32_t val) {
    for (int offset = 16; offset > 0; offset /= 2)
        val += __shfl_down_sync(0xFFFFFFFF, val, offset);

    return val;
}

__device__ inline uint16_t get_mask_for_range(int start_range, int end_range, int word_start_bit) {
    int word_end_bit = word_start_bit + 15;
    int isect_start = max(word_start_bit, start_range);
    int isect_end = min(word_end_bit, end_range);

    if (isect_start > isect_end)
        return 0;

    int len = isect_end - isect_start + 1;
    uint16_t partial_mask = (len == 16) ? 0xFFFF : ((1U << len) - 1);
    return partial_mask << (isect_start - word_start_bit);
}

__device__ inline uint16_t mask_range_512_cuda(const int offset, const int count) {
    const int lane_id = threadIdx.x & 31;
    uint16_t mask;

    if (count >= 512) {
        mask = 0xFFFF;
    } else {
        const int word_start_bit = lane_id * 16;
        int end = offset + count - 1;
        uint16_t m1;
        uint16_t m2 = 0;

        if (end < 512) {
            m1 = get_mask_for_range(offset, end, word_start_bit);
        } else {
            m1 = get_mask_for_range(offset, 511, word_start_bit);
            m2 = get_mask_for_range(0, end & 511, word_start_bit);
        }

        mask = m1 | m2;
    }

    __syncwarp();
    return mask;
}

__device__ inline uint16_t mask_range_512_cuda_rot(const int offset, const int count) {
    const int lane_id = threadIdx.x & 31;
    if (count >= 512)
        return 0xFFFF;

    int bits = count - (lane_id << 4);
    bits = max(0, min(bits, 16));
    const uint16_t mask = (bits == 16) ? 0xFFFF : ((1u << bits) - 1);
    const unsigned ofs = (0 - offset) & 511;
    const unsigned w = ofs >> 4;
    const unsigned b = ofs & 15;
    const uint16_t curr_word = __shfl_sync(0xFFFFFFFF, mask, (lane_id + w) & 31);
    const uint16_t next_word = __shfl_sync(0xFFFFFFFF, curr_word, (lane_id + 1) & 31);

    return (curr_word >> b) | (next_word << (16 - b));
}

__device__ inline uint16_t rotate_right_512(const uint16_t word, unsigned r) {
    const int lane_id = threadIdx.x & 31;
    r &= 511;
    const unsigned w = r >> 4;
    const unsigned b = r & 15;
    const uint16_t curr_word = __shfl_sync(0xFFFFFFFF, word, (lane_id + w) & 31);
    const uint16_t next_word = __shfl_sync(0xFFFFFFFF, curr_word, (lane_id + 1) & 31);

    return (curr_word >> b) | (next_word << (16 - b));
}

__device__ inline uint16_t rotate_right_512_deng(uint16_t word, unsigned r) {
    const int lane = threadIdx.x & 31;
    r &= 511;
    const unsigned w = r >> 4;
    const unsigned b = r & 15;
    uint16_t curr = __shfl_sync(0xffffffff, word, (lane + w) & 31);
    uint16_t next = __shfl_sync(0xffffffff, word, (lane + w + 1) & 31);

    return (curr >> b) | (next << (16 - b));
}

__device__ inline void print_512bit(uint16_t *a) {
    uint64_t *p = (uint64_t *)a;
    for (int i = 0; i < 8; i++)
        printf("%016llx ", p[i]);
    printf("\n");
}

__device__ inline void print_512float(__half *a) {
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 16; j++)
            printf("%f ", __half2float(a[i * 16 + j]));
        printf("\n");
    }
}

__device__ inline void print_512cn_msg(cn_msg *a) {
    printf("print min1_val\n");
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 16; j++)
            printf("%f ", __half2float(a[i * 16 + j].min1_val));
        printf("\n");
    }

    printf("print min2_val\n");
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 16; j++)
            printf("%f ", __half2float(a[i * 16 + j].min2_val));
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

__device__ inline bool codeword_bit_active(int col_idx, int bit_idx) {
    const int parity_start_col = decoder_global_info.cols - decoder_global_info.rows;
    if (col_idx < parity_start_col)
        return (col_idx * decoder_global_info.bits + bit_idx) < (decoder_global_info.bytes_of_userdata * 8);

    return ((col_idx - parity_start_col) * decoder_global_info.bits + bit_idx) < (decoder_global_info.bytes_of_parity * 8);
}

__device__ void inject_normal_errors(decoder_input_cw *input_cw) {
    const int lane_id = threadIdx.x & 31;
    const float rber = decoder_global_info.error_region_prob[0];
    const int state_idx = threadIdx.x + blockIdx.x * blockDim.x;
    curandState localState = d_states[state_idx];
    const auto idx_limit = decoder_global_info.cols * codeword::WORD_COUNT;

    for (int i = lane_id; i < idx_limit; i += WARP_SIZE) {
        const auto col_idx = i / codeword::WORD_COUNT;
        const auto word_idx = i % codeword::WORD_COUNT;
        uint64_t error_pattern = 0;
        uint64_t soft0_pattern = 0;
        uint64_t soft1_pattern = 0;

        for (unsigned bit_idx = 0; bit_idx < 64; bit_idx++) {
            const int bit_in_col = word_idx * 64 + bit_idx;
            if (!codeword_bit_active(col_idx, bit_in_col))
                continue;

            float vn = curand_uniform(&localState);
            if (vn < rber)
                error_pattern |= codec_support::UNIT << bit_idx;

            if (decoder_global_info.soft_bits == 1) {
                if ((decoder_global_info.error_region_prob[2] < vn) && (vn < decoder_global_info.error_region_prob[1])) {
                    soft0_pattern |= codec_support::UNIT << bit_idx;
                    soft1_pattern |= codec_support::UNIT << bit_idx;
                }
            } else if (decoder_global_info.soft_bits == 2) {
                if ((decoder_global_info.error_region_prob[2] < vn) && (vn < decoder_global_info.error_region_prob[1])) {
                    soft0_pattern |= codec_support::UNIT << bit_idx;
                    soft1_pattern |= codec_support::UNIT << bit_idx;
                } else if ((decoder_global_info.error_region_prob[4] < vn) && (vn < decoder_global_info.error_region_prob[3])) {
                    soft1_pattern |= codec_support::UNIT << bit_idx;
                } else if ((decoder_global_info.error_region_prob[6] < vn) && (vn < decoder_global_info.error_region_prob[5])) {
                    soft0_pattern |= codec_support::UNIT << bit_idx;
                }
            }
        }
        input_cw->hard.cols[col_idx][word_idx] = error_pattern;
        input_cw->soft0.cols[col_idx][word_idx] = soft0_pattern;
        input_cw->soft1.cols[col_idx][word_idx] = soft1_pattern;
    }

    d_states[state_idx] = localState;

    __syncwarp();

    return;
}

__device__ void inject_normal_errors2(decoder_input_cw *input_cw) {
    const int lane_id = threadIdx.x & 31;

    const int warp_id = threadIdx.x / 32;
    const int cw_id = blockIdx.x * NUM_WARPS + warp_id;

    const int state_idx = threadIdx.x + blockIdx.x * blockDim.x;
    curandState localState = d_states[state_idx];
    const auto idx_limit = decoder_global_info.cols * codeword::WORD_COUNT;

    __syncwarp();

    for (int i = lane_id; i < idx_limit; i += WARP_SIZE) {
        const auto col_idx = i / codeword::WORD_COUNT;
        const auto word_idx = i % codeword::WORD_COUNT;
        uint64_t error_pattern = 0;
        uint64_t soft0_pattern = 0;
        uint64_t soft1_pattern = 0;

        uint64_t bit = 1ULL;
        for (unsigned bit_idx = 0; bit_idx < 64; bit_idx += 2) {
            float2 r = curand_normal2(&localState);
            float vn0 = 1.0f + decoder_global_info.awgn_sigma * r.x;
            float vn1 = 1.0f + decoder_global_info.awgn_sigma * r.y;

#pragma unroll
            for (int k = 0; k < 2; k++) {
                const int bit_in_col = word_idx * 64 + bit_idx + k;
                if (!codeword_bit_active(col_idx, bit_in_col)) {
                    bit <<= 1;
                    continue;
                }

                float vn = (k == 0) ? vn0 : vn1;
                float avn = fabsf(vn);
                bool e = (vn < decoder_global_info.vref[0]);
                bool s1 = (vn < decoder_global_info.vref[3]) & (vn > decoder_global_info.vref[4]);
                bool c1 = (avn > decoder_global_info.vref[3]) & (avn < decoder_global_info.vref[5]);
                bool s0 = (avn < decoder_global_info.vref[1]) | c1;

                error_pattern |= (uint64_t)e * bit;
                soft1_pattern |= (uint64_t)s1 * bit;
                soft0_pattern |= (uint64_t)s0 * bit;
                bit <<= 1;
            }
        }

        input_cw->hard.cols[col_idx][word_idx] = error_pattern;
        input_cw->soft0.cols[col_idx][word_idx] = soft0_pattern;
        input_cw->soft1.cols[col_idx][word_idx] = soft1_pattern;
    }

    d_states[state_idx] = localState;

    return;
}

__device__ inline void update_node(int layer_pre, int col, decoder_workspace *cw_workspace, OptimizedSharedMemory *cw_shared_mem, int layer) {
    int sign_tmp;
    const int warp_id = threadIdx.x >> 5;
    const int lane_id = threadIdx.x & 31;
    const int cw_id = blockIdx.x * NUM_WARPS + warp_id;

    uint16_t mask_bit = 1;
    const char occupied = decoder_global_info.occupied[layer_pre][col];
    const char fade = decoder_global_info.fade[layer_pre][col];
    const bool last_row = (layer_pre == (decoder_global_info.rows - 1));
    const bool occu_last_row = occupied && last_row;
    const bool need_mask = last_row || fade;
    const int elem = decoder_global_info.element[layer_pre][col];
    __half cn_r_new_pre, cn_app_pre;
    int bit_idx, word_idx, bit_in_word;

    if (need_mask) {
        bit_idx = (lane_id + elem) & 511;
        word_idx = bit_idx >> 4;
        bit_in_word = bit_idx & 15;
    }

    for (int i = lane_id; i < P_SIZE; i += WARP_SIZE) {
        if (need_mask) {
            uint16_t word = decoder_global_info.mask[col][word_idx];
            mask_bit = (word >> bit_in_word) & 1;
            word_idx = (word_idx + 2) & 31;
        }

        __half cn_q_sel_pre = cw_workspace[cw_id].vn.cn_q_mem[col][i];
#ifdef USE_CN_MSG
        const cn_msg cn_c_sel_pre = cw_workspace[cw_id].vn.cn_c_mem[layer_pre][i];
#else
        cn_msg cn_c_sel_pre;
        cn_c_sel_pre.min1_val = cw_workspace[cw_id].vn.cn_c_mem_min1_val[layer_pre][i];
        cn_c_sel_pre.min2_val = cw_workspace[cw_id].vn.cn_c_mem_min2_val[layer_pre][i];
        uint8_t cn_c_sel_pre_min1_pos = cw_workspace[cw_id].vn.cn_c_mem_min1_pos[layer_pre][i];
        cn_c_sel_pre.sign_tot = (cn_c_sel_pre_min1_pos >> 7) ? 1 : -1;
        cn_c_sel_pre.min1_pos = cn_c_sel_pre_min1_pos & 127;
#endif

        sign_tmp = __half2float(cn_q_sel_pre) >= 0 ? 1 : -1;
        if (cn_c_sel_pre.min1_pos == col)
            cn_r_new_pre = cn_c_sel_pre.min2_val * __half(cn_c_sel_pre.sign_tot * sign_tmp);
        else
            cn_r_new_pre = cn_c_sel_pre.min1_val * __half(cn_c_sel_pre.sign_tot * sign_tmp);

        const bool cond1 = (occupied && !last_row) || (occu_last_row && mask_bit) || (fade && (!mask_bit));

        if (cond1)
            cn_app_pre = cn_r_new_pre + cn_q_sel_pre;
        else
            cn_app_pre = cn_q_sel_pre;

#ifdef FINITE_MODE
        cn_app_pre = Sat_Quan(cn_app_pre, decoder_global_info.finite_q_max, decoder_global_info.finite_q_min, decoder_global_info.finite_q_num, decoder_global_info.finite_f_num);
#endif

        cw_shared_mem->cn_app_pre[i] = cn_app_pre;
    }
}

__device__ inline void update_hd(int col, int layer, int layer_pre, short &hd_updated, decoder_workspace *cw_workspace, OptimizedSharedMemory *cw_shared_mem, int &shift_val1) {
    const int warp_id = threadIdx.x / 32;
    const int lane_id = threadIdx.x & 31;
    const int cw_id = blockIdx.x * NUM_WARPS + warp_id;
    const bool occu_cond = decoder_global_info.occupied[layer][col] && (layer == (decoder_global_info.rows - 1));
    const bool fade = decoder_global_info.fade[layer][col];
    int shift_val2;

    if (cw_shared_mem->dec_init[col] == 1) {
        shift_val1 = decoder_global_info.element[layer][col];
        shift_val2 = 0;
        if (lane_id == 0)
            cw_shared_mem->dec_init[col] = 0;
    } else {
        shift_val1 = (P_SIZE - decoder_global_info.element[layer_pre][col] + decoder_global_info.element[layer][col]) & P_SIZEm1;
        shift_val2 = P_SIZE - decoder_global_info.element[layer_pre][col];
    }

    uint16_t vn_dec_hd = 0;
    int base = lane_id * 16 + shift_val2;
#pragma unroll
    for (int i = 0; i < 16; i++) {
        int idx = base + i;
        if (idx >= P_SIZE)
            idx -= P_SIZE;

        unsigned short raw = *reinterpret_cast<unsigned short *>(&cw_shared_mem->cn_app_pre[idx]);

        vn_dec_hd |= ((raw >> 15)) << i;
    }

    cw_shared_mem->vn_dec_hd[lane_id] = vn_dec_hd;
    uint64_t *p = (uint64_t *)cw_shared_mem->vn_dec_hd;

    if (hd_updated == 0) {
        uint64_t diff = 0;
        if (lane_id < 8)
            diff = p[lane_id] ^ cw_workspace[cw_id].vn.dec_do_blk.cols[col][lane_id];
        unsigned mask = __ballot_sync(0xffffffff, (lane_id < 8) && (diff != 0));
        hd_updated = (mask != 0);
    }

    if (lane_id < 8)
        cw_workspace[cw_id].vn.dec_do_blk.cols[col][lane_id] = p[lane_id];

    vn_dec_hd = rotate_right_512_deng(vn_dec_hd, decoder_global_info.element[layer][col]);
    uint16_t rr = rotate_right_512_deng(decoder_global_info.mask[col][lane_id], decoder_global_info.element[layer][col]);
    if (occu_cond)
        vn_dec_hd &= rr;
    if (fade)
        vn_dec_hd &= ~rr;

    cw_shared_mem->layer_synd[lane_id] ^= vn_dec_hd;
    __syncwarp();
}

__device__ inline void update_node_next(int cir_cnt, int layer, int col, OptimizedSharedMemory *cw_shared_mem, decoder_workspace *cw_workspace, int shift_val1) {
    int sign_tmp;
    __half val_tmp;
    const int warp_id = threadIdx.x / 32;
    const int lane_id = threadIdx.x & 31;
    const int cw_id = blockIdx.x * NUM_WARPS + warp_id;
    const int elem = decoder_global_info.element[layer][col];
    uint16_t mask_bit = 1;
    const char occupied = decoder_global_info.occupied[layer][col];
    const char fade = decoder_global_info.fade[layer][col];
    const bool need_mask = (layer == (decoder_global_info.rows - 1)) || fade;
    const bool occu_cond = occupied && (layer == (decoder_global_info.rows - 1));
    __half cn_r_old, cn_q_updt_cur;
    uint32_t cn_q_sign;

    if (!need_mask) {
        for (int i = lane_id; i < P_SIZE; i += WARP_SIZE) {
            cn_msg cn_c_updt_cur = cw_workspace[cw_id].cn_c_updt_cur[i];
#ifdef USE_CN_MSG
            cn_msg cn_c_sel_cur = cw_workspace[cw_id].vn.cn_c_mem[layer][i];
#else
            cn_msg cn_c_sel_cur;
            cn_c_sel_cur.min1_val = cw_workspace[cw_id].vn.cn_c_mem_min1_val[layer][i];
            cn_c_sel_cur.min2_val = cw_workspace[cw_id].vn.cn_c_mem_min2_val[layer][i];
            uint8_t cn_c_sel_cur_min1_pos = cw_workspace[cw_id].vn.cn_c_mem_min1_pos[layer][i];
            cn_c_sel_cur.sign_tot = (cn_c_sel_cur_min1_pos >> 7) ? 1 : -1;
            cn_c_sel_cur.min1_pos = cn_c_sel_cur_min1_pos & 127;
#endif

            if (lane_id == 0)
                cn_q_sign = cw_workspace[cw_id].cn_q_sign[cir_cnt][i >> 5];
            cn_q_sign = __shfl_sync(0xffffffff, cn_q_sign, 0);

            int bit = (cn_q_sign >> lane_id) & 1;
            int sign = ((bit << 1) - 1) * cn_c_sel_cur.sign_tot;
            int use_min2 = (cn_c_sel_cur.min1_pos == col);

            __half diff = __hsub(cn_c_sel_cur.min2_val, cn_c_sel_cur.min1_val);
            __half val = cn_c_sel_cur.min1_val + __int2half_rn(use_min2) * diff;

            cn_r_old = __hmul(val, __int2half_rn(sign));
            const __half cn_app_cur = cw_shared_mem->cn_app_pre[(i + shift_val1) & P_SIZEm1];
            cn_q_updt_cur = cn_app_cur - cn_r_old;

#ifdef FINITE_MODE
            cn_q_updt_cur = Sat_Quan(cn_q_updt_cur, decoder_global_info.finite_q_max, decoder_global_info.finite_q_min, decoder_global_info.finite_q_num, decoder_global_info.finite_f_num);
#endif

            sign_tmp = (cn_q_updt_cur >= __half(0.0f)) ? 1 : -1;
            val_tmp = cn_q_updt_cur * __half(sign_tmp);
            cn_c_updt_cur.sign_tot *= sign_tmp;

            if (val_tmp < cn_c_updt_cur.min1_val) {
                cn_c_updt_cur.min2_val = cn_c_updt_cur.min1_val;
                cn_c_updt_cur.min1_val = val_tmp;
                cn_c_updt_cur.min1_pos = col;
            } else if (val_tmp < cn_c_updt_cur.min2_val) {
                cn_c_updt_cur.min2_val = val_tmp;
            }

            if (sign_tmp == -1)
                sign_tmp = 0;

            uint32_t warp_mask = __ballot_sync(0xffffffff, sign_tmp);
            if (lane_id == 0)
                cw_workspace[cw_id].cn_q_sign[cir_cnt][i >> 5] = warp_mask;

            cw_workspace[cw_id].vn.cn_q_mem[col][i] = cn_q_updt_cur;
            cw_workspace[cw_id].cn_c_updt_cur[i] = cn_c_updt_cur;
        }
    } else {
        for (int i = lane_id; i < P_SIZE; i += WARP_SIZE) {
            uint16_t bit_idx_in512 = (i + elem) & 511;
            uint16_t word_idx = bit_idx_in512 >> 4;
            uint16_t bit_idx_in_word = bit_idx_in512 & 15;
            mask_bit = (decoder_global_info.mask[col][word_idx] >> bit_idx_in_word) & 0x1;

            cn_msg cn_c_updt_cur = cw_workspace[cw_id].cn_c_updt_cur[i];
#ifdef USE_CN_MSG
            cn_msg cn_c_sel_cur = cw_workspace[cw_id].vn.cn_c_mem[layer][i];
#else
            cn_msg cn_c_sel_cur;
            cn_c_sel_cur.min1_val = cw_workspace[cw_id].vn.cn_c_mem_min1_val[layer][i];
            cn_c_sel_cur.min2_val = cw_workspace[cw_id].vn.cn_c_mem_min2_val[layer][i];
            uint8_t cn_c_sel_cur_min1_pos = cw_workspace[cw_id].vn.cn_c_mem_min1_pos[layer][i];
            cn_c_sel_cur.sign_tot = (cn_c_sel_cur_min1_pos >> 7) ? 1 : -1;
            cn_c_sel_cur.min1_pos = cn_c_sel_cur_min1_pos & 127;
#endif

            if (lane_id == 0)
                cn_q_sign = cw_workspace[cw_id].cn_q_sign[cir_cnt][i >> 5];
            cn_q_sign = __shfl_sync(0xffffffff, cn_q_sign, 0);

            int bit = (cn_q_sign >> lane_id) & 1;
            int sign = ((bit << 1) - 1) * cn_c_sel_cur.sign_tot;
            int use_min2 = (cn_c_sel_cur.min1_pos == col);
            __half diff = __hsub(cn_c_sel_cur.min2_val, cn_c_sel_cur.min1_val);
            __half val = cn_c_sel_cur.min1_val + __int2half_rn(use_min2) * diff;

            cn_r_old = __hmul(val, __int2half_rn(sign));
            const bool spec_cond = (occu_cond && !mask_bit) || (fade && mask_bit);
            const __half cn_app_cur = cw_shared_mem->cn_app_pre[(i + shift_val1) & P_SIZEm1];

            if (spec_cond)
                cn_r_old = 0;
            cn_q_updt_cur = cn_app_cur - cn_r_old;

#ifdef FINITE_MODE
            cn_q_updt_cur = Sat_Quan(cn_q_updt_cur, decoder_global_info.finite_q_max, decoder_global_info.finite_q_min, decoder_global_info.finite_q_num, decoder_global_info.finite_f_num);
#endif

            if (spec_cond) {
                sign_tmp = 1;
                val_tmp = 10000.0;
            } else {
                sign_tmp = (cn_q_updt_cur >= __half(0.0f)) ? 1 : -1;
                val_tmp = cn_q_updt_cur * __half(sign_tmp);
            }

            cn_c_updt_cur.sign_tot *= sign_tmp;

            if (val_tmp < cn_c_updt_cur.min1_val) {
                cn_c_updt_cur.min2_val = cn_c_updt_cur.min1_val;
                cn_c_updt_cur.min1_val = val_tmp;
                cn_c_updt_cur.min1_pos = col;
            } else if (val_tmp < cn_c_updt_cur.min2_val) {
                cn_c_updt_cur.min2_val = val_tmp;
            }

            if (sign_tmp == -1)
                sign_tmp = 0;

            uint32_t warp_mask = __ballot_sync(0xffffffff, sign_tmp);
            if (lane_id == 0)
                cw_workspace[cw_id].cn_q_sign[cir_cnt][i >> 5] = warp_mask;

            cw_workspace[cw_id].vn.cn_q_mem[col][i] = cn_q_updt_cur;
            cw_workspace[cw_id].cn_c_updt_cur[i] = cn_c_updt_cur;
        }
    }
}
