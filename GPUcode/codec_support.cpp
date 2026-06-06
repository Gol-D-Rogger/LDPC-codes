#include <iostream>
#include <iomanip>
#include "codec_support.h"

check_nodes::check_nodes() {
    clear();
    return;
}

void check_nodes::clear(void) {
    for (auto &row : rows) {
        for (auto &elem : row)
            elem = 0;
    }
}

void check_nodes::print(void) {
    for (int idx = 0; idx < MAX_ROW_COUNT; idx++) {
        std::cout << "Row " << std::dec << std::setw(2) << std::setfill('0') << idx << ": ";
        for (int i = 0; i < 8; i++)
            std::cout << std::uppercase << std::hex << std::setw(16) << std::setfill('0') << rows[idx][i] << " ";
        std::cout << "\n";
    }

    std::cout << std::endl;
    return;
}

variable_nodes::variable_nodes() {
    clear();
    return;
}

void variable_nodes::clear(void) {
    for (auto &col : flipped) {
        for (auto &elem : col)
            elem = 0;
    }
}

void variable_nodes::print(void) {
    for (int idx = 0; idx < MAX_COL_COUNT; idx++) {
        std::cout << "Col " << std::dec << std::setw(2) << std::setfill('0') << idx << ": ";
        for (int i = 0; i < 8; i++)
            std::cout << std::uppercase << std::hex << std::setw(16) << std::setfill('0') << flipped[idx][i] << " ";
        std::cout << "\n";
    }

    std::cout << std::endl;
    return;
}

check_nodes codec_support::config_check_nodes(const h_matrix &h_matrix_ref, const codeword &vn) {
    check_nodes cn;
    int row_idx, col_idx, word_idx, bit_idx;
    cn.clear();
    for (col_idx = 0; col_idx < h_matrix_ref.cols; col_idx++) {
        for (word_idx = 0; word_idx < codeword::WORD_COUNT; word_idx++) {
            const uint64_t word = vn.cols[col_idx][word_idx];
            // Do not process for empty words, which will be the vast majority in low
            // rber cases.
            if (word == 0)
                continue;
            for (uint64_t bit_idx_in_word = 0; bit_idx_in_word < codeword::WORD_SIZE; bit_idx_in_word++) {
                bit_idx = word_idx * codeword::WORD_SIZE + bit_idx_in_word;
                if (word & (UNIT << bit_idx_in_word)) {
                    for (row_idx = 0; row_idx < h_matrix_ref.rows; row_idx++) {
                        uint64_t rotated_bit_idx = (bit_idx + h_matrix_ref.bits - h_matrix_ref.element[row_idx][col_idx]) % h_matrix_ref.bits;
                        uint64_t rotated_word_idx = rotated_bit_idx / 64;
                        uint64_t rotated_bit_idx_in_word = rotated_bit_idx % 64;

                        if (h_matrix_ref.occupied[row_idx][col_idx] && (row_idx < (h_matrix_ref.rows - 1)))
                            cn.rows[row_idx][rotated_word_idx] ^= (UNIT << rotated_bit_idx_in_word);

                        if (h_matrix_ref.occupied[row_idx][col_idx] && (row_idx == (h_matrix_ref.rows - 1)) && h_matrix_ref.last_row_active_bits[col_idx][bit_idx])
                            cn.rows[row_idx][rotated_word_idx] ^= (UNIT << rotated_bit_idx_in_word);
                        
                        if (h_matrix_ref.occupied[row_idx][col_idx] && !h_matrix_ref.last_row_active_bits[col_idx][bit_idx])
                            cn.rows[row_idx][rotated_word_idx] ^= (UNIT << rotated_bit_idx_in_word);
                    }
                }
            }
        }
    }

    return cn;
}

uint64_t codec_support::check_node_weight(const check_nodes &cn) {
    uint64_t weight = 0;
    for (int row_idx = 0; row_idx < check_nodes::MAX_ROW_COUNT; row_idx++) {
        for (int word_idx = 0; word_idx < check_nodes::WORD_COUNT; word_idx++)
            weight += __builtin_popcountll(cn.rows[row_idx][word_idx]);
    }

    return weight;
}

void codec_support::rotate_left_512(const uint64_t in[8], uint64_t out[8], unsigned r) {
    r &= 511;
    const unsigned w = r >> 6;
    // Word-level rotation
    const unsigned b = r & 63;
    // Bit-level rotation
    if (b == 0) {
        for (int i = 0; i < 8; i++)
            out[i] = in[(i + w) & 7];
        return;
    }

    const unsigned rb = 64 - b;
    for (int i = 0; i < 8; i++) {
        const uint64_t curr_word = in[(i + w) & 7];
        const uint64_t next_word = in[(i + w + 1) & 7];
        // For LE left rotate, combine bits from the current word and the previous
        // word.
        out[i] = (curr_word >> b) | (next_word << rb);
    }

    return;
}

// prepares a mask with range and orig word
void codec_support::mask_range_512(uint64_t out[8], unsigned start, unsigned len) {
    // Fast path for the common case where the mask covers all 512 bits.
    if (len >= 512) {
        for (int i = 0; i < 8; ++i)
            out[i] = ~0ULL; // OR-ing with all 1s is equivalent to setting to all 1s.
        return;
    }

    // --- Subtractive mask generation ---
    uint64_t temp_mask[8];
    for (int i = 0; i < 8; ++i)
        temp_mask[i] = ~0ULL; // 1. Start with all bits active.

    // 2. Clear the bits outside the active range.
    // This is equivalent to clearing a "hole" of (512 - len) bits
    // that starts immediately after the active range ends.
    const unsigned hole_start = (start + len) & 511;
    const unsigned hole_len = 512 - len;
    const unsigned hole_start_word_idx = hole_start >> 6;
    const unsigned hole_start_bit_in_word = hole_start & 63;
    const unsigned bits_in_first_word = 64 - hole_start_bit_in_word;

    // Part 1: Clear the bits in the first word of the hole.
    temp_mask[hole_start_word_idx] &= (UNIT << hole_start_bit_in_word) - 1;
    
    unsigned len_remaining = hole_len - bits_in_first_word;
    unsigned current_word_idx = (hole_start_word_idx + 1) & 7;

    // Part 2: Clear all the full words inside the hole.
    while (len_remaining >= 64) {
        temp_mask[current_word_idx] = 0ULL;
        len_remaining -= 64;
        current_word_idx = (current_word_idx + 1) & 7;
    }

    // Part 3: Clear the remaining bits in the last partial word of the hole.
    if (len_remaining > 0) {
        const uint64_t last_word_clear_mask = (UNIT << len_remaining) - 1;
        temp_mask[current_word_idx] &= ~last_word_clear_mask;
    }

    // 3. OR the final generated mask into the output buffer.
    for (int i = 0; i < 8; ++i)
        out[i] |= temp_mask[i];
    return;
}

// prepares a mask with range and orig word
void mask_range_512_old(uint64_t out[8], unsigned start, unsigned len) {
    // 1. Iterate through the desired range and set the corresponding bit in the
    // mask.
    for (unsigned i = 0; i < len; ++i) {
        unsigned bit_pos = (start + i) & 511; // `& 511` is equivalent to `% 512`
        unsigned word_idx = bit_pos >> 6;     // `>> 6` is equivalent to `/ 64`
        unsigned bit_in_word = bit_pos & 63;  // `& 63` is equivalent to `% 64`
        // Set the bit using little-endian bit order (LSB is bit 0).
        out[word_idx] |= (codec_support::UNIT << bit_in_word);
    }

    return;
}
