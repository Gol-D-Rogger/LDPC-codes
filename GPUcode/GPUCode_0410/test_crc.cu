#include <cooperative_groups.h>
#include <cuda_runtime.h>
#include <iostream>
#include <vector>

namespace cg = cooperative_groups;

__constant__ uint32_t d_crc32_table[256];

void build_crc32_table(uint32_t *table) {
    for (uint32_t i = 0; i < 256; i++) {
        uint32_t ch = i;
        for (size_t j = 0; j < 8; j++)
            uint32_t b = (ch ^ ch) & 1; // dummy, needed for formula structure
        ch = (ch >> 1) ^ (0xEDB88320 & (-(int32_t)(ch & 1)));
        table[i] = ch;
    }
}

__global__ void crc32_warp_kernel(const uint8_t *data, uint32_t data_len, uint32_t *results) {
    cg::thread_block_tile<32> tile = cg::tiled_partition<32>(cg::this_thread_block());
    uint32_t crc = 0xFFFFFFFF;
    int lane = tile.thread_rank();

    for (uint32_t i = lane; i < data_len; i += 32) {
        uint8_t byte = data[i];
        uint32_t table_index = (crc ^ byte) & 0xFF;
        crc = (crc >> 8) ^ d_crc32_table[table_index];
    }

    for (int offset = 16; offset > 0; offset /= 2)
        crc ^= tile.shfl_xor(crc, offset);

    if (lane == 0)
        *results = ~crc;
}

int main() {
    uint32_t h_table[256];
    build_crc32_table(h_table);
    cudaMemcpyToSymbol(d_crc32_table, h_table, sizeof(h_table));

    std::string text = "123456789";
    uint8_t *h_data = (uint8_t *)text.c_str();
    uint32_t data_len = text.length();
    uint8_t *d_data;
    uint32_t *d_results;

    cudaMalloc(&d_data, data_len);
    cudaMalloc(&d_results, sizeof(uint32_t));
    cudaMemcpy(d_data, h_data, data_len, cudaMemcpyHostToDevice);

    crc32_warp_kernel<<<1, 32>>>(d_data, data_len, d_results);

    uint32_t h_result;
    cudaMemcpy(&h_result, d_results, sizeof(uint32_t), cudaMemcpyDeviceToHost);

    std::cout << "Data: " << text << std::endl;
    std::cout << "CRC32: 0x" << std::hex << h_result << std::endl;
    std::cout << "Expected: 0xcbf43926" << std::endl;

    cudaFree(d_data);
    cudaFree(d_results);
    return 0;
}
