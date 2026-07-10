#pragma once
#include <vector>
#include <array>
#include <cstdint>

struct codeword {
    static constexpr int WORD_COUNT = 8;
    // 8 words x 64 bits = 512
    static constexpr int WORD_SIZE = 64;
    static constexpr int MAX_COL_COUNT = 80;
    uint64_t cols[MAX_COL_COUNT][WORD_COUNT];
    void clear(void);
    void print(void);
};

struct decoder_input_cw {
    codeword hard;
    codeword soft0, soft1;
    decoder_input_cw();
    void clear(void);
};
